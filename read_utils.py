"""
read_utils.py — core per-read processing functions for SMF genome-wide analysis.

Extracted from genome_aligns.py and metaplot_around_bed.py to allow unit testing
without a live BAM file. Import these in scripts rather than duplicating logic.
"""

import numpy as np


def build_tally(xm_tag, read_positions, win_start, max_index):
    """Convert a Bismark XM tag and reference positions to a methylation tally.

    Parameters
    ----------
    xm_tag : str
        Bismark XM string (one character per aligned base).
        Alphabetic = CpG/CHG/CHH context; uppercase = unmethylated, lowercase = methylated.
        '.' or 'X' (non-context) are ignored.
    read_positions : array-like of int
        Reference genome positions for each base (pysam get_reference_positions()).
    win_start : int
        Genomic start of the analysis window; positions become relative to this.
    max_index : int
        Length of the output array (= window * 2 + 1).

    Returns
    -------
    numpy.ndarray of float, shape (max_index,), or None if no valid positions.
        NaN = no coverage; 0.0 = unmethylated context; 1.0 = methylated context.
    """
    xm_arr = np.asarray(list(xm_tag))
    rel_pos = np.asarray(read_positions, dtype=int) - win_start
    valid = (rel_pos >= 0) & (rel_pos < max_index) & np.char.isalpha(xm_arr)
    if not np.any(valid):
        return None
    tally = np.full(max_index, np.nan, dtype=float)
    tally[rel_pos[valid]] = 0.0
    tally[rel_pos[valid & np.char.islower(xm_arr)]] = 1.0
    return tally


def classify_read_strand(motif_strand, ys_tag):
    """Classify a read as same- or opposite-strand relative to a motif.

    Bismark non-directional PE assigns one of four YS strands:
      OT (original top) and CTOT (complementary to OT) → top-strand reads
      OB (original bottom) and CTOB (complementary to OB) → bottom-strand reads

    Parameters
    ----------
    motif_strand : str
        '+' or '-'
    ys_tag : str
        Bismark YS tag value: 'OT', 'CTOT', 'CTOB', or 'OB'.

    Returns
    -------
    'same' or 'opp'

    Raises
    ------
    ValueError for unrecognised inputs.
    """
    top = {'OT', 'CTOT'}
    bot = {'CTOB', 'OB'}
    if motif_strand == '+':
        if ys_tag in top: return 'same'
        if ys_tag in bot: return 'opp'
    elif motif_strand == '-':
        if ys_tag in top: return 'opp'
        if ys_tag in bot: return 'same'
    raise ValueError(f"Unrecognised motif_strand={motif_strand!r} or ys_tag={ys_tag!r}")


def annotate_positions(aligned_pairs, query_sequence, max_read_len=132):
    """Build a per-position annotation array from aligned pairs and query sequence.

    Separates bisulfite conversions (C→T, G→A) from true mismatches.
    Suitable for passing pysam's get_aligned_pairs(matches_only=True, with_seq=True).

    Parameters
    ----------
    aligned_pairs : iterable of (read_pos, ref_pos, ref_base)
        ref_base may be upper- or lower-case (pysam MD-tag convention).
    query_sequence : str
        Read sequence string (read.query_sequence).
    max_read_len : int
        Output array length.

    Returns
    -------
    numpy.ndarray dtype='U1', shape (max_read_len,)
        '.' = matched, non-informative
        'C' = C ref, C query  (unconverted cytosine)
        'T' = C ref, T query  (bisulfite-converted cytosine)
        'G' = G ref, G query  (unconverted guanine)
        'A' = G ref, A query  (converted guanine, bottom-strand reads)
        'M' = other mismatch
    """
    arr = np.full(max_read_len, '.', dtype='U1')
    for read_pos, _, ref_base in aligned_pairs:
        qbase = query_sequence[read_pos].upper()
        rbase = ref_base.upper()
        if rbase != qbase:
            arr[read_pos] = 'M'
        if rbase == 'C':
            if qbase == 'C':   arr[read_pos] = 'C'
            elif qbase == 'T': arr[read_pos] = 'T'
        if rbase == 'G':
            if qbase == 'G':   arr[read_pos] = 'G'
            elif qbase == 'A': arr[read_pos] = 'A'
    return arr


def is_chimeric(arr, threshold=2):
    """Return True if a read has both C→T and G→A conversions above threshold.

    Reads with abundant T *and* A annotations likely span OT+OB signal,
    indicating a chimeric PCR artifact.

    Parameters
    ----------
    arr : numpy.ndarray of dtype 'U1'
        Output of annotate_positions().
    threshold : int
        Minimum count of each conversion type to flag as chimeric.
    """
    return bool((np.sum(arr == 'T') > threshold) and (np.sum(arr == 'A') > threshold))


def minimal_unique_prefix_length(seqs, max_cycles):
    """Shortest prefix that uniquely distinguishes each index sequence.

    Parameters
    ----------
    seqs : list of str
        Already-deduplicated index sequences.
    max_cycles : int
        Maximum index read cycles (caps the search).

    Returns
    -------
    int : Minimum prefix length for unambiguous demultiplexing.
    """
    max_len = min(max_cycles, max(len(s) for s in seqs))
    for length in range(1, max_len + 1):
        if len({s[:length] for s in seqs}) == len(seqs):
            return length
    return max_len


# ---------------------------------------------------------------------------
# MD-based conversion calls (minibwa --meth / Bismark), used by the metaplots
# ---------------------------------------------------------------------------

def conversion_strand(read):
    """Which converted strand a read's bases come from, in BAM (reference) orientation.

    'OT': conversions show as ref C -> read T.  'OB': ref G -> read A.

    Bismark's XG tag (genome conversion, CT or GA) is used when present; it is
    set identically on both mates and also covers CTOT/CTOB in non-directional
    runs. Otherwise (minibwa --meth, directional) strand comes from the FLAG:
    single-end reads and R1 are OB when reverse; R2 is the opposite, because R2
    is sequenced from the complement of the same original strand.
    FLAG 99/147 -> OT, 83/163 -> OB, 0 -> OT, 16 -> OB.
    """
    if read.has_tag('XG'):
        return 'OT' if read.get_tag('XG') == 'CT' else 'OB'
    if read.is_paired and read.is_read2:
        return 'OT' if read.is_reverse else 'OB'
    return 'OB' if read.is_reverse else 'OT'


def project_to_window(read, win_start, win_len):
    """Place a read's aligned bases on window coordinates [win_start, win_start + win_len).

    Soft clips and insertions are skipped (matches_only); deletions and
    uncovered positions are 'N' with quality 0. Reference bases come from the
    read's MD tag, so no FASTA is needed.

    Returns
    -------
    (bases, quals, refs) : arrays of shape (win_len,), dtypes 'S1', int, 'S1'
    """
    bases = np.full(win_len, b'N', dtype='S1')
    quals = np.zeros(win_len, dtype=int)
    refs = np.full(win_len, b'N', dtype='S1')
    pairs = read.get_aligned_pairs(matches_only=True, with_seq=True)
    if not pairs:
        return bases, quals, refs
    qpos, rpos, rbase = zip(*pairs)
    qpos = np.asarray(qpos, dtype=int)
    rel = np.asarray(rpos, dtype=int) - win_start
    keep = (rel >= 0) & (rel < win_len)
    seq = np.frombuffer(read.query_sequence.encode('ascii'), dtype='S1')
    bases[rel[keep]] = seq[qpos[keep]]
    q = read.query_qualities
    quals[rel[keep]] = np.asarray(q, dtype=int)[qpos[keep]] if q is not None else 0
    refs[rel[keep]] = np.frombuffer(''.join(rbase).upper().encode('ascii'), dtype='S1')[keep]
    return bases, quals, refs


def merge_mates(mate1, mate2=None):
    """Merge one molecule's mates in reference coordinates (amplicon-smf
    vinsfan368 fork's merge_mates rules).

    Each mate is a (bases, quals, refs) tuple from project_to_window.
      1. Positions no mate covers are 'N'.
      2. Covered by one mate with base quality 0 -> 'N'.
      3. Covered by one mate with quality > 0 -> that base.
      4. Overlap, different qualities -> the higher-quality base.
      5. Overlap, same quality, same base -> that base.
      6. Overlap, same quality, different bases -> 'N'.
    With mate2=None (single-end, or the mate lies outside the window), only
    rules 1-3 apply.

    Returns
    -------
    (bases, refs) : arrays of shape (win_len,), dtype 'S1'
    """
    s0, q0, r0 = mate1
    if mate2 is None:
        return np.where(q0 > 0, s0, b'N'), r0
    s1, q1, r1 = mate2
    merged = np.full(s0.shape, b'N', dtype='S1')
    merged[q0 > q1] = s0[q0 > q1]
    merged[q1 > q0] = s1[q1 > q0]
    same = (q0 == q1) & (s0 == s1) & (q0 > 0)
    merged[same] = s0[same]
    refs = np.where(r0 != b'N', r0, r1)
    return merged, refs


def conversion_calls(strand, refs, bases):
    """Informative and converted positions for one molecule.

    OT: candidates are ref C, read C (retained) or T (converted).
    OB: candidates are ref G, read G (retained) or A (converted).
    Any other read base at a candidate (a mismatch, or 'N') is not informative.

    Returns
    -------
    (covered, converted) : boolean arrays
    """
    if strand == 'OB':
        candidate = refs == b'G'
        converted = candidate & (bases == b'A')
        retained = candidate & (bases == b'G')
    else:
        candidate = refs == b'C'
        converted = candidate & (bases == b'T')
        retained = candidate & (bases == b'C')
    return converted | retained, converted
