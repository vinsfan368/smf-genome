#!/usr/bin/env bash
set -euo pipefail

module load biology
module load bcl2fastq

############################
# CONSTANTS (EDIT THESE)   #
############################
MIN_BASES=8   # minimum index prefix length to enforce

# Where Illumina run folders live (you pass the *folder name* as an argument).
ILLUMINA_BASEDIR="$OAK/seqruns"   # e.g. "/oak/stanford/groups/.../seqruns"

###########
# PATHS   #
###########
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &>/dev/null && pwd )"

# Always look in script dir for samplesheet, and always write outputs there
SHEETDIR="${SCRIPT_DIR}"
OUTDIR="${SCRIPT_DIR}"

# Illumina data search location (run folders only)
SEARCHDIR="${ILLUMINA_BASEDIR}"

#############################
# USER ARG: RUN FOLDER NAME #
#############################
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <RUNFOLDER_NAME>  (e.g. $0 251126_M01234_0123_ABCDEFGH)" >&2
  exit 1
fi
RUNFOLDER_NAME="$1"

DATECODE="${RUNFOLDER_NAME%%_*}"

echo "Run folder name:     ${RUNFOLDER_NAME}"
echo "Datecode:            ${DATECODE}"
echo "Run search dir:      ${SEARCHDIR}"
echo "Samplesheet dir:     ${SHEETDIR}"
echo "Output directory:    ${OUTDIR}"

########################
# FIND RUN FOLDER
########################
INDIR="${SEARCHDIR}/${RUNFOLDER_NAME}"
if [[ ! -d "${INDIR}" ]]; then
  echo "ERROR: Run folder not found: ${INDIR}" >&2
  exit 1
fi

########################
# FIND SAMPLESHEET
########################
shopt -s nullglob
CANDIDATES=( "${SHEETDIR}/${DATECODE}"*.csv )
shopt -u nullglob

if (( ${#CANDIDATES[@]} == 0 )); then
  echo "ERROR: No samplesheet candidates found in ${SHEETDIR} matching ${DATECODE}*.csv" >&2
  exit 1
fi

SHEETS=()
for f in "${CANDIDATES[@]}"; do
  bn="$(basename "$f")"
  if [[ "${bn,,}" == *samplesheet* ]]; then
    SHEETS+=( "$f" )
  fi
done

if (( ${#SHEETS[@]} == 1 )); then
  SAMPLESHEET="${SHEETS[0]}"
elif (( ${#SHEETS[@]} > 1 )); then
  echo "ERROR: Multiple samplesheet candidates found in ${SHEETDIR}:" >&2
  printf '  %s\n' "${SHEETS[@]}" >&2
  exit 1
else
  if (( ${#CANDIDATES[@]} == 1 )); then
    SAMPLESHEET="${CANDIDATES[0]}"
  else
    echo "ERROR: Multiple ${DATECODE}*.csv candidates found, none contain 'samplesheet':" >&2
    printf '  %s\n' "${CANDIDATES[@]}" >&2
    exit 1
  fi
fi

echo "Using run folder:    ${INDIR}"
echo "Using samplesheet:   ${SAMPLESHEET}"

##################################
# PARSE READ CYCLES FROM [Reads] #
##################################
READ1CYCLES=$(awk -F '[,\t]' '$1=="Read1Cycles"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "${SAMPLESHEET}")
READ2CYCLES=$(awk -F '[,\t]' '$1=="Read2Cycles"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "${SAMPLESHEET}")
INDEX1CYCLES=$(awk -F '[,\t]' '$1=="Index1Cycles"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "${SAMPLESHEET}")
INDEX2CYCLES=$(awk -F '[,\t]' '$1=="Index2Cycles"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "${SAMPLESHEET}")

if [[ -z "${READ1CYCLES}" || -z "${READ2CYCLES}" || -z "${INDEX1CYCLES}" || -z "${INDEX2CYCLES}" ]]; then
  echo "ERROR: Failed to parse Read/Index cycles from [Reads] section" >&2
  exit 1
fi

echo "Read1Cycles  = ${READ1CYCLES}"
echo "Read2Cycles  = ${READ2CYCLES}"
echo "Index1Cycles = ${INDEX1CYCLES}"
echo "Index2Cycles = ${INDEX2CYCLES}"

#########################################
# PYTHON: MINIMAL UNIQUE INDEX PREFIXES #
#########################################
read MIN_I1 MIN_I2 <<< "$(
python3 - "${SAMPLESHEET}" "${INDEX1CYCLES}" "${INDEX2CYCLES}" "${MIN_BASES}" << 'PY'
import sys, re

ss_path = sys.argv[1]
index1_cycles = int(sys.argv[2])
index2_cycles = int(sys.argv[3])
MIN_BASES = int(sys.argv[4])

index1_seqs = []
index2_seqs = []

in_data = False
idx_idx1 = None
idx_idx2 = None

with open(ss_path, 'r', encoding='utf-8') as fh:
    for raw_line in fh:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith('[Data]'):
            in_data = True
            header_line = next(fh).strip()
            cols = re.split(r'[\t,]', header_line)
            norm = [c.strip().lower() for c in cols]
            try:
                idx_idx1 = norm.index('index')
            except ValueError:
                raise SystemExit("ERROR: Could not find 'Index' column in [Data] header")
            try:
                idx_idx2 = norm.index('index2')
            except ValueError:
                raise SystemExit("ERROR: Could not find 'index2'/'Index2' column in [Data] header")
            continue
        if not in_data:
            continue
        cols = re.split(r'[\t,]', line)
        if len(cols) <= max(idx_idx1, idx_idx2):
            continue
        i1 = cols[idx_idx1].strip()
        i2 = cols[idx_idx2].strip()
        if i1: index1_seqs.append(i1)
        if i2: index2_seqs.append(i2)

if not index1_seqs or not index2_seqs:
    raise SystemExit("ERROR: No index sequences found in [Data] section")

index1_unique = sorted(set(index1_seqs))
index2_unique = sorted(set(index2_seqs))

def minimal_unique_prefix_length(seqs, max_cycles):
    max_len = min(max_cycles, max(len(s) for s in seqs))
    for L in range(1, max_len + 1):
        if len({s[:L] for s in seqs}) == len(seqs):
            return L
    return max_len

min_i1 = max(MIN_BASES, min(minimal_unique_prefix_length(index1_unique, index1_cycles), index1_cycles))
min_i2 = max(MIN_BASES, min(minimal_unique_prefix_length(index2_unique, index2_cycles), index2_cycles))

print(min_i1, min_i2)
PY
)"

echo "Minimal (unique+>=${MIN_BASES}) prefix for Index  = ${MIN_I1} bases"
echo "Minimal (unique+>=${MIN_BASES}) prefix for index2 = ${MIN_I2} bases"

####################
# BUILD BASES MASK #
####################
if (( MIN_I1 > INDEX1CYCLES )); then
  echo "ERROR: MIN_I1 (${MIN_I1}) > Index1Cycles (${INDEX1CYCLES})" >&2
  exit 1
fi
if (( MIN_I2 > INDEX2CYCLES )); then
  echo "ERROR: MIN_I2 (${MIN_I2}) > Index2Cycles (${INDEX2CYCLES})" >&2
  exit 1
fi

I1_IGNORED=$(( INDEX1CYCLES - MIN_I1 ))
I2_IGNORED=$(( INDEX2CYCLES - MIN_I2 ))
I1_MASK="I${MIN_I1}$( (( I1_IGNORED > 0 )) && echo "n${I1_IGNORED}" || true )"
I2_MASK="I${MIN_I2}$( (( I2_IGNORED > 0 )) && echo "n${I2_IGNORED}" || true )"

BASEMASK="Y${READ1CYCLES},${I1_MASK},${I2_MASK},Y${READ2CYCLES}"
echo "Using bases mask: ${BASEMASK}"

##############################
# MAKE TRUNCATED SAMPLESHEET #
##############################
TRUNC_SAMPLESHEET="${OUTDIR}/${DATECODE}_samplesheet_trunc_i${MIN_I1}_${MIN_I2}.csv"
echo "Writing truncated samplesheet: ${TRUNC_SAMPLESHEET}"

python3 - "${SAMPLESHEET}" "${TRUNC_SAMPLESHEET}" "${MIN_I1}" "${MIN_I2}" << 'PY'
import sys, re

src = sys.argv[1]
dst = sys.argv[2]
min_i1 = int(sys.argv[3])
min_i2 = int(sys.argv[4])

with open(src, 'r', encoding='utf-8') as fin, open(dst, 'w', encoding='utf-8') as fout:
    in_data = False
    idx_idx1 = None
    idx_idx2 = None
    delim = ','

    lines = fin.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.rstrip('\n')

        if not in_data:
            fout.write(line)
            if stripped.startswith('[Data]'):
                in_data = True
                i += 1
                header_line = lines[i].rstrip('\n')
                fout.write(lines[i])
                delim = ',' if ',' in header_line else '\t'
                cols = [c.strip() for c in re.split(r'[,\t]', header_line)]
                norm = [c.lower() for c in cols]
                try:
                    idx_idx1 = norm.index('index')
                    idx_idx2 = norm.index('index2')
                except ValueError:
                    raise SystemExit("ERROR: Could not find Index/index2 in [Data] header")
        else:
            if not stripped or stripped.startswith('['):
                fout.write(line)
                in_data = False
                i += 1
                continue
            parts = stripped.split(delim)
            if len(parts) > max(idx_idx1, idx_idx2):
                if parts[idx_idx1]: parts[idx_idx1] = parts[idx_idx1][:min_i1]
                if parts[idx_idx2]: parts[idx_idx2] = parts[idx_idx2][:min_i2]
            fout.write(delim.join(parts) + '\n')
        i += 1
PY

#################
# RUN bcl2fastq #
#################
mkdir -p "${OUTDIR}/fastq"

bcl2fastq \
  --create-fastq-for-index-reads \
  --barcode-mismatches 0 \
  --use-bases-mask "${BASEMASK}" \
  --runfolder-dir "${INDIR}" \
  --output-dir "${OUTDIR}/fastq" \
  --sample-sheet "${TRUNC_SAMPLESHEET}" \
  --no-lane-splitting

########################
# WRITE DOWNSTREAM samplesheet.tsv
########################
SAMPLESHEET_TSV="${OUTDIR}/${DATECODE}_samplesheet.tsv"
echo "Writing downstream samplesheet: ${SAMPLESHEET_TSV}"

python3 - "${SAMPLESHEET}" "${SAMPLESHEET_TSV}" "${OUTDIR}" "${DATECODE}" << 'PY'
import sys, re, os

orig_ss = sys.argv[1]
out_tsv = sys.argv[2]
outdir  = sys.argv[3]
datecode = sys.argv[4]

fastq_dir = os.path.join(outdir, "fastq")

amplicon_fa = "/oak/stanford/groups/wjg/bgrd/projects/smf/241023_P092_deaminase/opJS45.amplicon.fa"
experiment = "deaminase"
filter_contigs = "FALSE"
no_endog_meth = "TRUE"
ignore_bounds = "FALSE"

rows = []
in_data = False
idx_sample_id = None
idx_sample_name = None
delim = ','

with open(orig_ss, 'r', encoding='utf-8') as fh:
    lines = fh.readlines()

i = 0
row_index = 0
while i < len(lines):
    line = lines[i]
    stripped = line.rstrip('\n')

    if not in_data:
        if stripped.startswith('[Data]'):
            in_data = True
            i += 1
            header_line = lines[i].rstrip('\n')
            delim = ',' if ',' in header_line else '\t'
            cols = [c.strip() for c in re.split(r'[,\t]', header_line)]
            norm = [c.lower() for c in cols]
            try:
                idx_sample_id = norm.index("sample_id")
            except ValueError:
                raise SystemExit("ERROR: Could not find Sample_ID in [Data] header")
            try:
                idx_sample_name = norm.index("sample_name")
            except ValueError:
                idx_sample_name = idx_sample_id
        i += 1
        continue

    if not stripped or stripped.startswith('['):
        break

    parts = [p.strip() for p in stripped.split(delim)]
    if len(parts) <= max(idx_sample_id, idx_sample_name):
        i += 1
        continue

    row_index += 1
    sample_id = parts[idx_sample_id]
    sample_name = parts[idx_sample_name]
    if not sample_id:
        i += 1
        continue

    rows.append((
        sample_name,
        os.path.join(fastq_dir, f"{sample_id}_S{row_index}_R1_001.fastq.gz"),
        os.path.join(fastq_dir, f"{sample_id}_S{row_index}_R2_001.fastq.gz"),
        amplicon_fa, experiment, filter_contigs, no_endog_meth, ignore_bounds,
    ))
    i += 1

with open(out_tsv, 'w', encoding='utf-8') as out:
    out.write("\t".join(["sample_name","fastq_R1","fastq_R2","amplicon_fa",
                         "experiment","filter_contigs","no_endog_meth","ignore_bounds"]) + "\n")
    for r in rows:
        out.write("\t".join(r) + "\n")
PY
