# smf-genome

Genome-wide single-molecule footprinting (SMF) pipeline for Bismark non-directional paired-end WGBS data.

## Pipeline

```
[0] demux               bcl2fastq with auto-computed bases mask
 ↓
[1] trim                TrimGalore (stringency=3, 9 bp 5'-clip, paired FastQC)
 ↓
[2] align               Bismark BT2 PE (--non_directional --strandID)
 ↓
[3] dedup               deduplicate_bismark
 ↓
[4] sort / index        samtools sort + index → ready for pysam random access
```

## Setup

```bash
mamba env create -f environment.yaml
conda activate smf-genome
```

## Usage

### Demultiplex (step 0)

Place the Illumina samplesheet (`<DATECODE>*samplesheet*.csv`) in the repo directory, then:

```bash
./smf_pipeline.sh demux 251215_NB551514_0395_AHHKH7BGYX
```

This calls `scripts/demux_mask_idx_bases.sh`, auto-detects the minimal unique index
prefix, builds the bases mask, and writes FASTQs to `fastq/`.

### Snakemake pipeline (steps 1–4)

1. Edit `samples.tsv` — one row per sample with absolute paths to R1/R2 FASTQs.
2. Edit `config.yaml` — set `bismark_index` and `dedup_bin` for your system.
3. Run:

```bash
# Dry run
snakemake -n

# Local (16 cores)
snakemake -c16

# SLURM
snakemake --executor cluster-generic \
  --cluster-generic-submit-cmd "sbatch --cpus-per-task={threads} --mem=32G -t 4:00:00"
```

To rerun from a specific step (e.g. after changing dedup settings):

```bash
snakemake --forcerun dedup -c16
```

### Manual step-by-step

```bash
./smf_pipeline.sh run sample_R1.fastq.gz sample_R2.fastq.gz [outdir]
```

## Analysis

### Deamination metaplots

`scripts/metaplot_around_bed.py` — per-position deamination fraction around
the centres of BED-like regions, from coordinate-sorted, indexed minibwa
`--meth` or Bismark BAMs (MD tags required).

```bash
python scripts/metaplot_around_bed.py aligned/*.bam \
    --regions K562-CTCF-ENCFF119XFJ-vs-FIMO-intersection2.gz --out-dir metaplots
# zoomed in, with a sequence logo from the reference underneath
python scripts/metaplot_around_bed.py aligned/*.bam \
    --regions K562-CTCF-ENCFF119XFJ-vs-FIMO-intersection2.gz --out-dir metaplots \
    --half-width 20 --logo
```

| Option | Description |
|---|---|
| `--half-width` | bp either side of each region centre to scan and plot (default 500) |
| `--logo` | add a sequence-logo subplot (one FASTA pass over every region) |
| `--ylim LO HI` | fix the y-axis range (default: auto-scale) |
| `--fasta` | reference FASTA for `--logo` (default: hg38-male) |
| `--recompute` | ignore the cached `.npz` |

- Plots are in motif orientation: `-` regions are reversed (and their logo
  reverse-complemented), so x runs 5'→3' along each region's strand.
- Mates are merged per molecule (`read_utils.merge_mates`, amplicon-smf
  rules), so bases covered by both mates count once.
- Converted strand: Bismark `XG` tag if present, else FLAG + R1/R2
  (FLAG 99/147 → C→T, 83/163 → G→A).
- Region layouts (4/6/7/11/22 columns) are documented in `regions.py`.

### Sequence logos

```bash
python sequence_logos.py REGIONS [...] --out-dir logos [--by-strand] [--half-window 10]
python scripts/plot_strand_split_logos.py REGIONS [...] --out-dir logos [--half-window 15]
```

`plot_strand_split_logos.py` draws `+` / `-` / pooled logos side by side and
reports the best `+`/`-` lag, to catch BEDs whose centre is off by a uniform
offset. `motif_logo.py` has the subplot helper the metaplot uses.

### Other

`scripts/genome_aligns.py` — plots methylation profiles around CTCF motifs
and per-read conversion patterns. Set `BAM` and `regions` paths at the top,
then run directly.

Core per-read logic lives in `read_utils.py`, region parsing in `regions.py`
(importable, tested).

## Tests

```bash
pytest
```

## Config reference

| Key | Description |
|---|---|
| `samples` | Path to samples TSV (sample, R1, R2) |
| `bismark_index` | Bismark BT2 genome index directory |
| `dedup_bin` | Full path to `deduplicate_bismark` binary |
| `trim_cores` | Cores for TrimGalore |
| `bismark_parallel` | `--parallel` value for Bismark (total threads ≈ ×4) |
