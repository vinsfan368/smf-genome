#!/usr/bin/env python3
"""
Deamination metaplots around BED-like regions, from minibwa --meth or Bismark
alignments. Per-base conversion is decoded from each read's own MD tag, so no
genome FASTA is needed for the signal itself (only for the optional logo).

Molecules, not reads: within each region, a proper pair's two mates are merged
in reference coordinates with read_utils.merge_mates (amplicon-smf fork rules,
quality-based), so bases covered by both mates are counted once. Converted
strand per molecule comes from read_utils.conversion_strand (Bismark XG tag if
present, else FLAG + R1/R2).

Read filters: primary, not QC-fail, MAPQ >= MIN_MAPQ; for paired reads also
proper pair and mate MAPQ >= MIN_MAPQ (MQ tag, when the aligner writes one).

Plots are in motif orientation: '-' regions are reversed (and their logo
reverse-complemented) so x runs 5'->3' along each region's own strand, x<0
upstream. Region layouts are described in regions.py; region centre is
(start + end) // 2. Each BAM x region file x half-width scan is cached to .npz.

Usage:
  python scripts/metaplot_around_bed.py BAM [BAM ...] --regions REGIONS --out-dir DIR \\
      [--half-width 500] [--logo] [--ylim LO HI] [--recompute]
"""
import argparse
import os
import random
import sys

import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from read_utils import conversion_strand, project_to_window, merge_mates, conversion_calls
from regions import region_layout, parse_region_line, iter_region_tuples, open_regions
import sequence_logos as sl
import motif_logo as ml


MIN_MAPQ = 20
MAX_REGIONS = 5_000_000
STRANDS = ('OT', 'OB')
CACHE_VERSION = 'merged-mates-motif-v1'  # bump when compute_signal's output changes meaning


def _label(motif_strand):
    return 'pos' if motif_strand == '+' else 'neg'


def _passes_filters(read):
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
        return False
    if read.mapping_quality < MIN_MAPQ:
        return False
    if read.is_paired:
        if not read.is_proper_pair:
            return False
        if read.has_tag('MQ') and read.get_tag('MQ') < MIN_MAPQ:
            return False
    return True


def compute_signal(bam_path, regions_path, half_width):
    """Per-region mean conversion by position (equal weight per region), split
    by motif strand x converted strand, in motif orientation ('-' regions
    reversed)."""
    if not os.path.exists(bam_path + '.bai'):
        pysam.index(bam_path)

    field_indices, analysis_label = region_layout(regions_path)
    win_len = 2 * half_width + 1

    conv_sum = {s: {r: np.zeros(win_len) for r in STRANDS} for s in ('+', '-')}
    conv_cnt = {s: {r: np.zeros(win_len) for r in STRANDS} for s in ('+', '-')}

    with open_regions(regions_path) as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.startswith('#')]
    if len(lines) > MAX_REGIONS:
        lines = random.sample(lines, MAX_REGIONS)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for line in lines:
            contig, start, end, motif_strand = parse_region_line(line, field_indices)
            win_start = (int(start) + int(end)) // 2 - half_width
            try:
                region = bam.fetch(contig, win_start, win_start + win_len)
            except ValueError as err:
                if 'invalid contig' in str(err) or 'start out of range' in str(err):
                    continue
                raise

            molecules = {}
            for read in region:
                if _passes_filters(read):
                    molecules.setdefault(read.query_name, []).append(read)

            reg_sum = {r: np.zeros(win_len) for r in STRANDS}
            reg_cnt = {r: np.zeros(win_len) for r in STRANDS}
            for reads in molecules.values():
                strands = {conversion_strand(r) for r in reads}
                if len(reads) > 2 or len(strands) != 1:
                    continue
                rs = strands.pop()
                bases, refs = merge_mates(*(project_to_window(r, win_start, win_len) for r in reads))
                covered, converted = conversion_calls(rs, refs, bases)
                reg_cnt[rs][covered] += 1.0
                reg_sum[rs][converted] += 1.0

            for rs in STRANDS:
                c = reg_cnt[rs]
                covered_bins = c > 0
                if not covered_bins.any():
                    continue
                reg_mean = np.where(covered_bins, reg_sum[rs] / np.where(covered_bins, c, 1.0), 0.0)
                if motif_strand == '-':  # windows are symmetric, so reversing maps x -> -x
                    reg_mean, covered_bins = reg_mean[::-1], covered_bins[::-1]
                conv_sum[motif_strand][rs] += reg_mean
                conv_cnt[motif_strand][rs] += covered_bins.astype(float)

    return {
        'x': np.arange(-half_width, half_width + 1),
        'analysis_label': analysis_label,
        **{f'conv_sum_{_label(ms)}_{rs}': conv_sum[ms][rs] for ms in ('+', '-') for rs in STRANDS},
        **{f'conv_cnt_{_label(ms)}_{rs}': conv_cnt[ms][rs] for ms in ('+', '-') for rs in STRANDS},
    }


def _check_pooling(name, a_sum, a_cnt, b_sum, b_cnt, label_a, label_b):
    """Warn when the two profiles being pooled strongly ANTI-correlate.

    Both are in motif orientation, so they should agree; strong
    anti-correlation means one is mirror-imaged (e.g. a BED whose strand
    column is wrong), and averaging them cancels the signal: the pooled page
    goes flat while each per-strand page looks clean. Near-zero correlation is
    normal on low-conversion data, so only r < -0.3 warns.
    Prints a warning; never raises."""
    both = (a_cnt > 0) & (b_cnt > 0)
    n = int(both.sum())
    if n < 20:
        return None
    fa = a_sum[both] / a_cnt[both]
    fb = b_sum[both] / b_cnt[both]
    if fa.std() < 1e-9 or fb.std() < 1e-9:
        return None
    r = float(np.corrcoef(fa, fb)[0, 1])
    if r < -0.3:
        print(f'WARNING: {name}: its two components ({label_a} vs {label_b}) '
              f'ANTI-correlate (r={r:+.2f}, n={n} shared bins). Averaging them '
              f'cancels the signal -- the pooled trace is not meaningful here; '
              f'read the per-strand pages.', file=sys.stderr)
    return r


def plot_signal(data, bam_base, out_dir, regions_path, ylim=None, logo=False, fasta_path=sl.GENOME_FASTA):
    """One PDF: pooled same/opposite-strand page, then one page per motif strand
    x converted strand. ylim: (lo, hi) to fix the y-axis (with dashed y=0/y=1
    lines), or None to auto-scale. logo: add a sequence-logo subplot from the
    reference (one FASTA pass over every region)."""
    x = data['x']
    half_width = int(x.max())
    analysis_label = str(data['analysis_label'])
    xlabel = f"position relative to {analysis_label} center (bp, 5'->3' along each region's strand)"
    conv_sum = {ms: {rs: data[f'conv_sum_{_label(ms)}_{rs}'] for rs in STRANDS} for ms in ('+', '-')}
    conv_cnt = {ms: {rs: data[f'conv_cnt_{_label(ms)}_{rs}'] for rs in STRANDS} for ms in ('+', '-')}

    seq_probs = None
    if logo:
        # compute_base_probs reverse-complements '-' regions; the pooled logo
        # is the n-weighted mean of the two strands' logos
        by_strand = {'+': [], '-': []}
        for region in iter_region_tuples(regions_path):
            if region[3] in by_strand:
                by_strand[region[3]].append(region)
        fasta = pysam.FastaFile(fasta_path)
        seq_probs, n_used = {}, {}
        for ms in ('+', '-'):
            seq_probs[ms], n_used[ms], _ = sl.compute_base_probs(iter(by_strand[ms]), fasta, half_width)
        n_tot = n_used['+'] + n_used['-']
        seq_probs['pooled'] = ((n_used['+'] * seq_probs['+'] + n_used['-'] * seq_probs['-'])
                               / n_tot if n_tot else seq_probs['+'])

    def _plot(ax, y, **kw):
        """Drop NaN bins rather than letting them break the line: each trace is
        only defined at its own candidate base, so gaps are interleaved everywhere."""
        finite = np.isfinite(y)
        ax.plot(x[finite], y[finite], **kw)

    def _page(pdf, traces, title, logo_key):
        if seq_probs is not None:
            fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
        else:
            fig, ax = plt.subplots()
        for y, kw in traces:
            _plot(ax, y, **kw)
        ax.axvline(0, color='k', linewidth=0.5, linestyle='--')
        ax.set_xlim(-half_width, half_width)
        if ylim is not None:
            ax.set_ylim(*ylim)
            ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
            ax.axhline(1, color='k', linewidth=0.5, linestyle='--')
        ax.set_ylabel('deamination fraction')
        ax.set_title(title)
        if any('label' in kw for _, kw in traces):
            ax.legend()
        if seq_probs is not None:
            ml.draw_logo_ax(ax2, x, seq_probs[logo_key], xlabel)
        else:
            ax.set_xlabel(xlabel)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    def _mean(s, c):
        with np.errstate(invalid='ignore'):
            return np.where(c > 0, s / c, np.nan)

    same = (('+', 'OT'), ('-', 'OB'))
    opp = (('+', 'OB'), ('-', 'OT'))
    for nm, (ka, kb) in (('same strand', same), ('opposite strand', opp)):
        _check_pooling(f'{bam_base} "{nm}"',
                       conv_sum[ka[0]][ka[1]], conv_cnt[ka[0]][ka[1]],
                       conv_sum[kb[0]][kb[1]], conv_cnt[kb[0]][kb[1]],
                       f'motif{ka[0]}/{ka[1]}', f'motif{kb[0]}/{kb[1]}')
    same_mean = _mean(sum(conv_sum[m][r] for m, r in same), sum(conv_cnt[m][r] for m, r in same))
    opp_mean = _mean(sum(conv_sum[m][r] for m, r in opp), sum(conv_cnt[m][r] for m, r in opp))

    os.makedirs(out_dir, exist_ok=True)
    pdf_path = os.path.join(out_dir, f'{bam_base}_{analysis_label}_{half_width}bp_signal.pdf')
    with PdfPages(pdf_path) as pdf:
        _page(pdf, [(same_mean, dict(linewidth=0.8, label='same strand', alpha=0.6)),
                    (opp_mean, dict(linewidth=0.8, label='opposite strand', alpha=0.6))],
              bam_base, 'pooled')
        for ms in ('+', '-'):
            for rs in STRANDS:
                if conv_cnt[ms][rs].any():
                    _page(pdf, [(_mean(conv_sum[ms][rs], conv_cnt[ms][rs]), {})],
                          f'{bam_base}\nmotif strand: {ms}   read is {rs}', ms)
    print(f'saved: {pdf_path}')


def _params_key(bam_path, regions_path, half_width):
    """Cache key for the .npz: includes the BAM's and the regions file's own
    mtime+size, so editing either in place invalidates the cache."""
    st = os.stat(bam_path)
    rst = os.stat(regions_path)
    return (f'{CACHE_VERSION}|{half_width}|{MIN_MAPQ}|{MAX_REGIONS}|{regions_path}'
            f'|{st.st_mtime}|{st.st_size}|{rst.st_mtime}|{rst.st_size}')


def process_bam(bam_path, regions_path, out_dir, half_width=500, recompute=False, ylim=None, logo=False,
                fasta_path=sl.GENOME_FASTA):
    _, analysis_label = region_layout(regions_path)
    bam_base = os.path.splitext(os.path.basename(bam_path))[0]
    cache_path = os.path.join(out_dir, f'{bam_base}_{analysis_label}_{half_width}bp_signal_data.npz')
    key = _params_key(bam_path, regions_path, half_width)

    data = None
    if not recompute and os.path.exists(cache_path):
        try:
            cached = dict(np.load(cache_path, allow_pickle=True))
            if str(cached.get('params_key')) == key:
                data = cached
        except Exception:
            data = None

    if data is None:
        data = compute_signal(bam_path, regions_path, half_width)
        data['params_key'] = np.array(key)
        os.makedirs(out_dir, exist_ok=True)
        np.savez(cache_path, **data)
        print(f'cached: {cache_path}')
    else:
        print(f'loaded cache: {cache_path}')

    plot_signal(data, bam_base, out_dir, regions_path, ylim=ylim, logo=logo, fasta_path=fasta_path)


def main(argv=None):
    p = argparse.ArgumentParser(description='Deamination metaplots around BED-like regions.')
    p.add_argument('bams', nargs='+', help='coordinate-sorted, indexed BAM file(s)')
    p.add_argument('--regions', required=True, help='region file; see regions.py for supported layouts')
    p.add_argument('--out-dir', required=True, help='output directory for PDFs and the .npz cache')
    p.add_argument('--half-width', type=int, default=500,
                   help='bp either side of each region centre to scan and plot (default: %(default)s)')
    p.add_argument('--logo', action='store_true',
                   help='add a sequence-logo subplot from the reference (readable only for small half-widths)')
    p.add_argument('--ylim', type=float, nargs=2, metavar=('LO', 'HI'),
                   help='fix the y-axis range (default: auto-scale)')
    p.add_argument('--fasta', default=sl.GENOME_FASTA,
                   help='reference FASTA, used only for --logo (default: %(default)s)')
    p.add_argument('--recompute', action='store_true', help='ignore cached .npz data')
    args = p.parse_args(argv)

    for bam_path in args.bams:
        try:
            process_bam(bam_path, args.regions, args.out_dir, half_width=args.half_width,
                        recompute=args.recompute, ylim=args.ylim, logo=args.logo, fasta_path=args.fasta)
        except Exception as exc:
            print(f'FAILED: {bam_path}: {exc}', file=sys.stderr)


if __name__ == '__main__':
    main()
