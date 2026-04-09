#!/usr/bin/env python3
import argparse
import pysam
import gzip
import os

import pysam

import numpy as np
from collections import Counter
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.gridspec as grd

from read_utils import build_tally, classify_read_strand, annotate_positions, is_chimeric

# {genome: {read: strand}} nested dict to store how reads are aligned
BAM = '/Users/vinsfan/Desktop/longer_reads/trimgalore_extra_clip/aligns/dedup/2J_bismark_bt2_pe.deduplicated.bam'
UNCONV_FA = '/Users/vinsfan/Desktop/_greenleaf_lab/hg38/hg38-male.fa'
MAX_READ_LEN = 132
MAX_N_READS = 100_000
MAX_N_MISMATCHES = 10


if __name__ == "__main__":
    out = os.path.splitext(BAM)[0] + '.sorted.bam'
    pysam.sort('-o', out, BAM)
    pysam.index(out)
    bam = pysam.AlignmentFile(out, 'rb')
    regions = '/Users/vinsfan/Desktop/_greenleaf_lab/hg38/K562-CTCF-ENCFF119XFJ-vs-FIMO-intersection2.gz'
    window = 500
    max_index = (window * 2) + 1
    conv_fracs = {'+': {'OT': [], "CTOT": [], "CTOB": [], "OB": []},
                  '-': {'OT': [], "CTOT": [], "CTOB": [], "OB": []}}
    with gzip.open(regions, 'rt') as fh:
        for line in fh:
            _, _, contig, start, end, motif_strand, _, _, _, _, _ = line.split('\t')
            center = (int(start) + int(end)) // 2
            win_start = center - window
            win_end = center + window
            try:
                region = bam.fetch(contig, win_start, win_end)
            except ValueError as err:
                if "invalid contig" in str(err):
                    continue
                if 'start out of range' in str(err):
                    continue
                else:
                    raise

            region_conv = {'+': {'OT': [], "CTOT": [], "CTOB": [], "OB": []},
                           '-': {'OT': [], "CTOT": [], "CTOB": [], "OB": []}}
            for read in region:
                read_pos = np.asarray(read.get_reference_positions(), dtype=int)
                xm_tag = read.get_tag('XM')
                tally = build_tally(xm_tag, read_pos, win_start, max_index)
                if tally is None:
                    continue
                region_conv[motif_strand][read.get_tag('YS')].append(tally)

            for ms, dict2 in region_conv.items():
                for rs, fracs in dict2.items():
                    if len(fracs) > 0:
                        conv_fracs[ms][rs].extend(fracs)

    # Aggregate same- vs opposite-strand profiles
    same_strand = []
    opp_strand = []
    for motif_strand, dict2 in conv_fracs.items():
        for read_strand, fracs_list in dict2.items():
            fracs = np.asarray(fracs_list)
            if motif_strand == '-':
                fracs = fracs[:, ::-1]
            dest = same_strand if classify_read_strand(motif_strand, read_strand) == 'same' else opp_strand
            dest.append(fracs)

    same_strand = np.concatenate(same_strand, axis=0)
    opp_strand = np.concatenate(opp_strand, axis=0)
    same_mean = np.nanmean(same_strand, axis=0)
    opp_mean = np.nanmean(opp_strand, axis=0)
    plt.plot(1-same_mean, linewidth=0.8, label='same strand')
    plt.plot(1-opp_mean, linewidth=0.8, label='opposite strand')
    plt.ylim((0.6, 1.02))
    plt.ylabel("1 - deamination")
    plt.legend()
    plt.savefig(f"signal_around_CTCF.png", dpi=600)
    plt.show()

    # Per-strand-type plots
    for motif_strand, dict2 in conv_fracs.items():
        for read_strand, fracs_list in dict2.items():
            arr = np.asarray(fracs_list, dtype=float)
            mean = np.nanmean(arr, axis=0)
            plt.plot(1 - mean)
            plt.title(f"motif strand: {motif_strand}\nread is {read_strand}")
            plt.ylim((0.25, 1.02))
            plt.savefig(f"{motif_strand}_{read_strand}.png", dpi=600)
            plt.show()
            plt.close()


# Iterate to collect methylation tags and template insert lengths
tags = defaultdict(list)
template_lengths = defaultdict(list)
with pysam.AlignmentFile(BAM, 'rb') as bam:
    for i, read in enumerate(bam):
        if i % 2 != 0:  # only look at first strand
            continue
        read_strand = read.get_tag('YS')
        tlen = read.template_length
        if read_strand == 'OT' or read_strand == 'CTOB':
            if tlen < 0:
                continue
        elif read_strand == 'CTOT' or read_strand == 'OB':
            if tlen > 0:
                continue
        template_lengths[read_strand].append(read.template_length)
        xm_tag = read.get_tag('XM')
        if len(xm_tag) == 132:
            tags[read_strand].append(list(xm_tag))

# Plot conversion fraction as a function of read position
for key in tags.keys():
    arr = np.asarray(tags[key])
    print(arr.shape)
    aligns, positions = arr.shape
    conv_arr = np.full_like(arr, fill_value=np.nan, dtype=float)
    conv_arr[np.char.isalpha(arr)] = 0
    conv_arr[np.char.islower(arr)] = 1
    mean = np.nanmean(conv_arr, axis=0)
    coverage = np.sum(~np.isnan(conv_arr), axis=0)

    fig = plt.figure()
    gs = grd.GridSpec(4, 1, figure=fig)
    axs = [fig.add_subplot(gs[:3]), fig.add_subplot(gs[-1])]
    axs[0].plot(mean)
    axs[0].set_ylim((0, 1.0))
    axs[0].set_ylabel("fraction converted")
    axs[1].bar(range(coverage.shape[0]), coverage)
    axs[1].hlines(y=aligns, xmin=0, xmax=positions, linestyles='--')
    axs[1].set_ylabel("informative positions")
    axs[1].set_xlabel("position along R1")
    axs[0].set_title(f"Conversion and coverage along aligned R1s\n{key} strand")
    plt.tight_layout()
    plt.savefig(f"{key}_strand_alignment_along.png", dpi=600)
    plt.show()

# Plot tlen hists
for key, val in template_lengths.items():
    plt.hist(val, bins=50)
    xlim = (0, 1000) if key in ('OT', 'CTOB') else (-1000, 0)
    plt.title(f"insert lengths\nR1 is {key} strand")
    plt.xlim(xlim)
    plt.savefig(f"tlen_{key}.png", dpi=600)
    plt.show()
    plt.close()


# Collect manually annotated mutations
arrs = []
query_strands = []
query_lens = []
with pysam.AlignmentFile(BAM, 'rb') as bam:
    for i, read in enumerate(bam):
        if i % 2 != 0:
            continue
        read_strand = read.get_tag('YS')
        tlen = read.template_length
        if read_strand == 'OT' or read_strand == 'CTOB':
            if tlen < 0:
                continue
        elif read_strand == 'CTOT' or read_strand == 'OB':
            if tlen > 0:
                continue
        arr = annotate_positions(
            read.get_aligned_pairs(matches_only=True, with_seq=True),
            read.query_sequence,
            max_read_len=MAX_READ_LEN,
        )
        if is_chimeric(arr):
            continue
        if np.sum(arr == 'M') <= MAX_N_MISMATCHES:
            query_strands.append(read_strand)
            query_lens.append(read.query_length)
            arrs.append(arr)
    arrs = np.stack(arrs)
    query_strands = np.stack(query_strands, dtype=str)
    query_lens = np.stack(query_lens, dtype=int)

    scores = np.full_like(arrs, np.nan, dtype=float)
    scores[arrs == 'C'] = 0
    scores[arrs == 'G'] = 0
    scores[arrs == 'T'] = 1
    scores[arrs == 'A'] = 1

    plt.plot(np.nanmean(scores, axis=0))
    plt.show()
