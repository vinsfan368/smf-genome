#!/usr/bin/env python3
"""
sequence_logos.py — raw-probability sequence logos around region centres.

Letter height at each position is the observed base fraction (heights sum to 1
per position), no background normalisation. Sequence comes from the reference
FASTA and is reverse-complemented for '-' regions, so logos reflect region
(motif) orientation, not genome orientation. Region files are read with
regions.py, the same parser the metaplots use.

Also the shared helpers behind the metaplot logo subplot (compute_base_probs,
letter glyphs); motif_logo.draw_logo_ax draws them onto an existing Axes.

Usage:
  python sequence_logos.py REGIONS [REGIONS ...] --out-dir DIR [--by-strand] [--half-window 10]
"""

import argparse
import os
import sys

import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D
from matplotlib.patches import PathPatch

from regions import iter_region_tuples

GENOME_FASTA = '/oak/stanford/groups/wjg/vinsfan/hg38_refs/hg38-male-bismark-indexes-BT2/hg38-male.fa'
BASES = 'ACGT'
BASE_IDX = {b: i for i, b in enumerate(BASES)}

_LETTER_COLORS = {'A': '#109648', 'C': '#255C99', 'G': '#F7B32B', 'T': '#D62839'}
_FONT = FontProperties(family='monospace', weight='bold')
_COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


def _revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _clean_basename(path):
    base = os.path.basename(path)
    for suffix in ('.bed.gz', '.100radius', '.bed'):
        if base.endswith(suffix):
            return base[:-len(suffix)]
    return base


def _idx_lookup():
    lut = np.full(256, -1, dtype=np.int8)
    for i, b in enumerate(BASES):
        lut[ord(b)] = i
    return lut


_IDX_LOOKUP = _idx_lookup()


def probs_from_counts(counts, totals):
    with np.errstate(invalid='ignore'):
        return np.where(totals[:, None] > 0, counts / totals[:, None], 0.0)


def compute_base_counts(region_iter, fasta, half_window):
    """Per-position base counts over +/-half_window around each region's centre.

    region_iter yields (chrom, start, end, strand); '-' regions are
    reverse-complemented. Returns (counts[win_len, 4], totals[win_len], n_used, n_skipped)."""
    win_len = 2 * half_window + 1
    counts = np.zeros((win_len, 4))
    totals = np.zeros(win_len)
    pos_range = np.arange(win_len)

    n_used = 0
    n_skipped = 0
    for chrom, start, end, strand in region_iter:
        center = (int(start) + int(end)) // 2
        win_start = center - half_window
        win_end = center + half_window + 1
        try:
            seq = fasta.fetch(chrom, win_start, win_end).upper()
        except (KeyError, ValueError):
            n_skipped += 1
            continue
        if len(seq) != win_len:
            n_skipped += 1
            continue
        if strand == '-':
            seq = _revcomp(seq)
        n_used += 1

        arr = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        base_idx = _IDX_LOOKUP[arr]
        valid = base_idx >= 0
        counts[pos_range[valid], base_idx[valid]] += 1
        totals[pos_range[valid]] += 1

    return counts, totals, n_used, n_skipped


def compute_base_probs(region_iter, fasta, half_window):
    """Returns (probs[2*half_window+1, 4], n_used, n_skipped)."""
    counts, totals, n_used, n_skipped = compute_base_counts(region_iter, fasta, half_window)
    return probs_from_counts(counts, totals), n_used, n_skipped


def compute_base_probs_by_strand(regions_path, fasta, half_window):
    """'+' and '-' regions in separate accumulators instead of pooled -- a
    forward/reverse handling bug then shows up as a frame shift between the two
    strands' logos instead of being smeared into one pooled plot.
    Returns ({strand: probs}, {strand: n_used}, n_skipped)."""
    by_strand = {'+': [], '-': []}
    n_skipped = 0
    for region in iter_region_tuples(regions_path):
        if region[3] in by_strand:
            by_strand[region[3]].append(region)
        else:
            n_skipped += 1
    probs, n_used = {}, {}
    for strand, regs in by_strand.items():
        probs[strand], n_used[strand], skipped = compute_base_probs(iter(regs), fasta, half_window)
        n_skipped += skipped
    return probs, n_used, n_skipped


def plot_probability_track(ax, x, probs):
    """4-line (A/C/G/T) base-probability-by-position track, for use as a
    subplot underneath a signal panel that shares the same x-axis."""
    for i, b in enumerate(BASES):
        ax.plot(x, probs[:, i], linewidth=0.8, label=b, color=_LETTER_COLORS[b])
    ax.set_ylim(0, 1)
    ax.set_ylabel('base probability')
    ax.legend(fontsize='small', ncol=4, loc='upper right')


def _normalized_letter_path(letter):
    """TextPath for `letter` remapped so its bounding box is exactly the unit
    square [0,1]x[0,1] -- lets every instance be placed with one scale+translate."""
    tp = TextPath((0, 0), letter, size=1, prop=_FONT)
    bbox = tp.get_extents()
    return tp.transformed(Affine2D().translate(-bbox.x0, -bbox.y0).scale(1 / bbox.width, 1 / bbox.height))


_LETTER_PATHS = {b: _normalized_letter_path(b) for b in BASES}


def _letter_patch(letter, x0, y0, width, height, color):
    transform = Affine2D().scale(width, height).translate(x0, y0)
    return PathPatch(_LETTER_PATHS[letter].transformed(transform), color=color, lw=0)


def draw_logo(probs, out_path, title, n_regions):
    half_window = (probs.shape[0] - 1) // 2
    fig, ax = plt.subplots(figsize=(8, 3))
    x = np.arange(-half_window, half_window + 1)
    for col, xi in enumerate(x):
        order = np.argsort(probs[col])  # ascending prob -> tallest letter drawn last, on top
        y0 = 0.0
        for idx in order:
            h = probs[col, idx]
            if h <= 0:
                continue
            letter = BASES[idx]
            ax.add_patch(_letter_patch(letter, xi - 0.5, y0, 1.0, h, _LETTER_COLORS[letter]))
            y0 += h
    ax.axvline(0, color='k', linewidth=0.5, linestyle='--')
    ax.set_xlim(-half_window - 0.5, half_window + 0.5)
    ax.set_ylim(0, 1.02)
    ax.set_xticks(x)
    ax.set_xlabel('position relative to region center (bp)')
    ax.set_ylabel('probability')
    ax.set_title(f'{title}\n(n={n_regions} regions)')
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description='Raw-probability sequence logos around region centres.')
    p.add_argument('regions', nargs='+', help='region file(s); see regions.py for supported layouts')
    p.add_argument('--out-dir', required=True, help='output directory for logo PDFs')
    p.add_argument('--fasta', default=GENOME_FASTA, help='reference FASTA (default: %(default)s)')
    p.add_argument('--half-window', type=int, default=10, help='bp either side of the centre (default: %(default)s)')
    p.add_argument('--by-strand', action='store_true', help="separate logos for '+' and '-' regions")
    args = p.parse_args(argv)

    fasta = pysam.FastaFile(args.fasta)
    os.makedirs(args.out_dir, exist_ok=True)
    hw = args.half_window
    for bed_path in args.regions:
        name = _clean_basename(bed_path)
        try:
            if args.by_strand:
                probs, n_used, n_skipped = compute_base_probs_by_strand(bed_path, fasta, hw)
            else:
                probs, n_used, n_skipped = compute_base_probs(iter_region_tuples(bed_path), fasta, hw)
        except Exception as exc:
            print(f'FAILED: {bed_path}: {exc}', file=sys.stderr)
            continue

        if args.by_strand:
            for strand, label in (('+', 'plus_strand'), ('-', 'minus_strand')):
                out_path = os.path.join(args.out_dir, f'{name}_{label}_center{hw}bp_probability_logo.pdf')
                draw_logo(probs[strand], out_path, f'{name} ({strand} strand)', n_used[strand])
                print(f'saved: {out_path}  (used {n_used[strand]}, skipped {n_skipped})')
        else:
            out_path = os.path.join(args.out_dir, f'{name}_center{hw}bp_probability_logo.pdf')
            draw_logo(probs, out_path, name, n_used)
            print(f'saved: {out_path}  (used {n_used}, skipped {n_skipped})')


if __name__ == '__main__':
    main()
