"""
read_utils.py — core per-read processing functions for SMF genome-wide analysis.

Extracted from genome_aligns.py to allow unit testing without a live BAM file.
Import these in genome_aligns.py rather than duplicating logic.
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
