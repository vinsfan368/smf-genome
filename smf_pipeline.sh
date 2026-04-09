#!/usr/bin/env bash
set -euo pipefail

# SMF sequencing pipeline
# Steps: [0] demux  [1] trim  [2] align  [3] dedup

############################
# CONSTANTS (EDIT THESE)   #
############################
BISMARK_IDX="${OAK:-/oak/stanford/groups/wjg}/${USER}/hg38_refs/hg38-male-bismark-indexes-BT2"
DEDUP_BIN="/Users/vinsfan/miniforge3/envs/smf_py3_v7/bin/deduplicate_bismark"
TRIM_CORES=4
BISMARK_PARALLEL=4

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &>/dev/null && pwd )"

############################
# USAGE                    #
############################
usage() {
  cat >&2 <<EOF
Usage: $0 <subcommand> [args...]

Subcommands:
  demux <RUNFOLDER_NAME>
      Demultiplex via bcl2fastq (calls scripts/demux_mask_idx_bases.sh).
      Writes FASTQs to ./fastq/ and a downstream samplesheet TSV.

  trim <R1.fastq.gz> <R2.fastq.gz> [outdir]
      TrimGalore: stringency=3, 9 bp 5' clip, paired, FastQC.
      Outputs to <outdir>/trimgalore/  (default: same dir as R1).

  align <R1_val_1.fq.gz> <R2_val_2.fq.gz> [outdir]
      Bismark PE, non-directional, --strandID, minins=76, maxins=1000.
      Outputs to <outdir>/aligns/.

  dedup <bismark_pe.bam>
      deduplicate_bismark on a Bismark PE BAM.

  run <R1.fastq.gz> <R2.fastq.gz> [outdir]
      Full pipeline: trim → align → dedup.
EOF
  exit 1
}

[[ $# -lt 1 ]] && usage
CMD="$1"; shift

############################
# SUBCOMMANDS              #
############################

cmd_demux() {
  [[ $# -lt 1 ]] && { echo "demux: missing RUNFOLDER_NAME" >&2; exit 1; }
  "${SCRIPT_DIR}/scripts/demux_mask_idx_bases.sh" "$1"
}

cmd_trim() {
  [[ $# -lt 2 ]] && { echo "trim: missing R1 R2" >&2; exit 1; }
  local r1; r1="$(realpath "$1")"
  local r2; r2="$(realpath "$2")"
  local outdir="${3:-$(dirname "$r1")}"
  mkdir -p "${outdir}/trimgalore"
  trim_galore \
    --fastqc \
    --stringency 3 \
    --cores "${TRIM_CORES}" \
    --clip_R1 9 \
    --clip_R2 9 \
    --paired \
    -o "${outdir}/trimgalore" \
    "${r1}" "${r2}"
}

cmd_align() {
  [[ $# -lt 2 ]] && { echo "align: missing trimmed R1 R2" >&2; exit 1; }
  local r1; r1="$(realpath "$1")"
  local r2; r2="$(realpath "$2")"
  local outdir="${3:-$(dirname "$r1")}"
  mkdir -p "${outdir}/aligns"
  bismark \
    "${BISMARK_IDX}" \
    -1 "${r1}" \
    -2 "${r2}" \
    --parallel "${BISMARK_PARALLEL}" \
    --non_directional \
    --minins 76 \
    --maxins 1000 \
    --strandID \
    -o "${outdir}/aligns"
}

cmd_dedup() {
  [[ $# -lt 1 ]] && { echo "dedup: missing BAM" >&2; exit 1; }
  "${DEDUP_BIN}" "$1"
}

cmd_run() {
  [[ $# -lt 2 ]] && { echo "run: missing R1 R2" >&2; exit 1; }
  local r1; r1="$(realpath "$1")"
  local r2; r2="$(realpath "$2")"
  local outdir="${3:-$(dirname "$r1")}"

  echo "[1/3] Trim..."
  cmd_trim "${r1}" "${r2}" "${outdir}"

  # TrimGalore names: {stem}_val_1.fq.gz / {stem}_val_2.fq.gz
  local r1_stem; r1_stem="$(basename "${r1%.fastq.gz}")"
  local trimmed_r1="${outdir}/trimgalore/${r1_stem}_val_1.fq.gz"
  local trimmed_r2="${outdir}/trimgalore/$(basename "${r2%.fastq.gz}")_val_2.fq.gz"

  echo "[2/3] Align..."
  cmd_align "${trimmed_r1}" "${trimmed_r2}" "${outdir}"

  local bam="${outdir}/aligns/${r1_stem}_val_1_bismark_bt2_pe.bam"

  echo "[3/3] Dedup..."
  cmd_dedup "${bam}"

  echo "Done: ${bam%.bam}.deduplicated.bam"
}

############################
# DISPATCH                 #
############################
case "${CMD}" in
  demux) cmd_demux "$@" ;;
  trim)  cmd_trim  "$@" ;;
  align) cmd_align "$@" ;;
  dedup) cmd_dedup "$@" ;;
  run)   cmd_run   "$@" ;;
  *)     echo "Unknown subcommand: ${CMD}" >&2; usage ;;
esac
