"""
motif_logo.py — letter-stack base-probability logo drawn as a subplot under a signal panel.

The track spans the same +/-half_window as the signal panel above it (not the
fixed +/-10bp of sequence_logos.py's standalone logos). Probabilities depend
only on the region file and the genome FASTA, not on any BAM, so get_probs can
cache them to .npz and reuse them across datasets.
"""

import os

import numpy as np
import pysam

from regions import iter_region_tuples
from sequence_logos import (
    GENOME_FASTA, BASES, _clean_basename, _LETTER_COLORS, _letter_patch, compute_base_probs,
)


def get_probs(bed_path, half_window, cache_dir, fasta_path=GENOME_FASTA):
    """Returns (x, probs[2*half_window+1, 4], n_regions) for bed_path, cached
    to cache_dir and invalidated on the BED file's own mtime."""
    cache_path = os.path.join(cache_dir, f'{_clean_basename(bed_path)}_center{half_window}bp_probs.npz')
    bed_mtime = os.stat(bed_path).st_mtime
    if os.path.exists(cache_path):
        cached = np.load(cache_path)
        if float(cached['bed_mtime']) == bed_mtime:
            return cached['x'], cached['probs'], int(cached['n_used'])

    fasta = pysam.FastaFile(fasta_path)
    probs, n_used, _ = compute_base_probs(iter_region_tuples(bed_path), fasta, half_window)
    x = np.arange(-half_window, half_window + 1)
    os.makedirs(cache_dir, exist_ok=True)
    np.savez(cache_path, x=x, probs=probs, n_used=n_used, bed_mtime=bed_mtime)
    return x, probs, n_used


def draw_logo_ax(ax, x, probs, xlabel):
    """Same letter-stack logo as sequence_logos.draw_logo, drawn onto an
    existing Axes so it can sit as a subplot under a signal panel."""
    for col, xi in enumerate(x):
        order = np.argsort(probs[col])  # ascending prob -> tallest letter on top
        y0 = 0.0
        for idx in order:
            h = probs[col, idx]
            if h <= 0:
                continue
            letter = BASES[idx]
            ax.add_patch(_letter_patch(letter, xi - 0.5, y0, 1.0, h, _LETTER_COLORS[letter]))
            y0 += h
    ax.axvline(0, color='k', linewidth=0.5, linestyle='--')
    ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('motif\nprobability')
