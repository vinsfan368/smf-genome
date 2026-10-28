#!/usr/bin/env python3
"""
Three-panel probability logo per region file: '+' regions, '-' regions, and the
two pooled.

Each per-strand logo can be crisp while the pooled one is blurred, because
pooling is the n-weighted mean of the two PWMs. If the BED's coordinates put
the motif off (start+end)//2 by a uniform amount in genomic space, the '-'
strand reverse complement negates that offset, so the two strands land on
opposite sides of x=0 and their average smears. The printed best +/- lag
quantifies the shift.

Regions are parsed with regions.py and oriented by
sequence_logos.compute_base_probs, so the panels show exactly what the
metaplot draws under its signal.

Usage:
  python scripts/plot_strand_split_logos.py REGIONS [REGIONS ...] --out-dir DIR [--half-window 15]
"""
import argparse
import os
import sys

import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sequence_logos as sl
from regions import iter_region_tuples, label_from_basename

INK, INK2, SURFACE = '#0b0b0b', '#52514e', '#fcfcfb'


def bits(p):
    return float(((p * np.log2(np.clip(p, 1e-9, 1))).sum(1) + 2).sum())


def strand_pwms(path, fa, hw):
    tup = {'+': [], '-': []}
    for region in iter_region_tuples(path):
        if region[3] in tup:
            tup[region[3]].append(region)
    out = {}
    for k, regions in (('plus', tup['+']), ('minus', tup['-']), ('pooled', tup['+'] + tup['-'])):
        probs, n, _ = sl.compute_base_probs(iter(regions), fa, hw)
        out[k] = (probs, n)
    return out


def draw(ax, probs, hw):
    x = np.arange(-hw, hw + 1)
    for col, xi in enumerate(x):
        y0 = 0.0
        for i in np.argsort(probs[col]):  # ascending -> tallest on top
            h = probs[col, i]
            if h <= 0:
                continue
            b = sl.BASES[i]
            ax.add_patch(sl._letter_patch(b, xi - 0.5, y0, 1.0, h, sl._LETTER_COLORS[b]))
            y0 += h
    ax.axvline(0, color=INK, lw=1.0, ls='--', zorder=5)
    ax.set_xlim(-hw - 0.5, hw + 0.5)
    ax.set_ylim(0, 1.02)
    ax.set_facecolor(SURFACE)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=14)
    ax.set_yticks([0, 0.5, 1])


def main(argv=None):
    p = argparse.ArgumentParser(description="'+' / '-' / pooled probability logos per region file.")
    p.add_argument('regions', nargs='+', help='region file(s); see regions.py for supported layouts')
    p.add_argument('--out-dir', required=True, help='output directory for PDF/PNG logos')
    p.add_argument('--fasta', default=sl.GENOME_FASTA, help='reference FASTA (default: %(default)s)')
    p.add_argument('--half-window', type=int, default=15, help='bp either side of the centre (default: %(default)s)')
    args = p.parse_args(argv)

    hw = args.half_window
    fa = pysam.FastaFile(args.fasta)
    os.makedirs(args.out_dir, exist_ok=True)
    for path in args.regions:
        name = label_from_basename(path)
        P = strand_pwms(path, fa, hw)
        fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True, facecolor=SURFACE)
        for ax, key, lab in zip(axes, ('plus', 'minus', 'pooled'), ('+ strand', '− strand', 'pooled')):
            probs, n = P[key]
            draw(ax, probs, hw)
            ax.set_title(f'{lab}      {bits(probs):.2f} bits      n={n:,}',
                         fontsize=17, color=INK, loc='left', pad=8)
        axes[-1].set_xlabel('position relative to motif centre (bp)', fontsize=17, color=INK2)
        axes[1].set_ylabel('probability', fontsize=17, color=INK2)
        fig.tight_layout()
        for ext in ('pdf', 'png'):
            fig.savefig(os.path.join(args.out_dir, f'{name}_strand_split_logo.{ext}'),
                        dpi=130, facecolor=SURFACE)
        plt.close(fig)

        # the quantitative companion to the picture
        pp, pm, po = P['plus'][0], P['minus'][0], P['pooled'][0]
        lags = {}
        for lag in range(-4, 5):
            a, b = ((pp[:2 * hw + 1 - lag], pm[lag:]) if lag >= 0
                    else (pp[-lag:], pm[:2 * hw + 1 + lag]))
            lags[lag] = (a * b).sum() / len(a)
        best = max(lags, key=lags.get)
        n_tot = P['plus'][1] + P['minus'][1]
        w = (P['plus'][1] * pp + P['minus'][1] * pm) / n_tot
        realigned = (P['plus'][1] * pp + P['minus'][1] * np.roll(pm, -best, axis=0)) / n_tot
        print(f'{name}: + {bits(pp):6.2f}  - {bits(pm):6.2f}  pooled {bits(po):6.2f} bits'
              f'   | best +/- lag = {best:+d}bp'
              f'   | pooled==mean? {np.abs(po - w).max() < 1e-12}'
              f'   | realigned pooled {bits(realigned):6.2f} bits')
    print(f'wrote *_strand_split_logo.pdf/.png in {args.out_dir}')


if __name__ == '__main__':
    main()
