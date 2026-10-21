#!/bin/bash

set -uo pipefail

# Input and Output
# Input: Skera BAMs, primers.fasta, whitelist
NAS_IN="/.../"
REF_IN="/.../"

NAS_OUT="/.../"
# LOCAL working directory (fast disk)
WORKDIR="/.../"
# Reference genome (local install)
REF_FA="/.../refdata-gex-GRCh38-2024-A/fasta/genome.fa"
# -----------------------------------------------------------------------------
PRIMERS="primers.fasta"
WL="GEMX_3p_whitelist.txt"
WL_RC="GEMX_3p_whitelist_RC.txt"
THREADS=16
SORT_THREADS=8

if [ $

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate isoseq_env || { echo "conda env isoseq_env not found"; exit 1; }

LOGDIR="${NAS_OUT}/logs_isoseq"
mkdir -p "$WORKDIR" "$NAS_OUT" "$LOGDIR"
cd "$WORKDIR" || { echo "Cannot enter $WORKDIR"; exit 1; }
FAILED="${NAS_OUT}/failed_isoseq.log"; : > "$FAILED"

rsync -ah "${REF_IN}/${PRIMERS}" . || { echo "Missing ${REF_IN}/${PRIMERS}"; exit 1; }
if [ -f "${REF_IN}/${WL_RC}" ]; then
  rsync -ah "${REF_IN}/${WL_RC}" .
else
  rsync -ah "${REF_IN}/${WL}" . || { echo "Missing ${REF_IN}/${WL}"; exit 1; }
  rev "$WL" | tr 'ACGT' 'TGCA' > "$WL_RC"
  rsync -ah --ignore-existing "$WL_RC" "${REF_IN}/"
fi

fetch() {
  local f=$1
  [ -f "$f" ] && return 0
  for d in "$NAS_OUT" "$NAS_IN"; do
    if [ -f "${d}/${f}" ]; then
      rsync -ah --ignore-missing-args "${d}/${f}" "${d}/${f}.pbi" ./ && return 0
    fi
  done
  return 1
}

step() {
  local lib=$1 name=$2 input=$3; shift 3
  local marker="${LOGDIR}/${lib}.${name}.done"
  if [ -f "$marker" ]; then echo "    [skip] $name"; return 0; fi
  fetch "$input" || { echo "    [FAILED] $name: input $input not found locally or on NAS"; return 1; }
  echo "    [$(date +%H:%M:%S)] $name"
  if "$@" > "${LOGDIR}/${lib}.${name}.log" 2>&1; then
    touch "$marker"
  else
    echo "    [FAILED] $name (see ${LOGDIR}/${lib}.${name}.log)"
    return 1
  fi
}

for LIB in "${LIBS[@]}"; do
  echo "=========================================="
  echo " Iso-Seq3 pipeline: $LIB"
  echo "=========================================="

  # Input Skera file
  IN=""
  for cand in "skera.${LIB}.bam" "skera.${LIB}_cDNA.bam"; do
    if [ -f "$cand" ] || [ -f "${NAS_IN}/${cand}" ]; then IN="$cand"; break; fi
  done
  if [ -z "$IN" ]; then
    echo "    no Skera file for $LIB in $NAS_IN, skipping"; echo "$LIB no_skera_input" >> "$FAILED"; continue
  fi
  P="skera.${LIB}"

  step "$LIB" lima "$IN" \
    lima --per-read --isoseq --log-level INFO -j $THREADS \
      "$IN" "$PRIMERS" "${P}.oriented.bam" \
    || { echo "$LIB lima" >> "$FAILED"; continue; }

  step "$LIB" tag "${P}.oriented.5p--3p.bam" \
    isoseq tag --design T-12U-16B --log-level INFO -j $THREADS \
      "${P}.oriented.5p--3p.bam" "${P}.tagged.bam" \
    || { echo "$LIB tag" >> "$FAILED"; continue; }

  step "$LIB" refine "${P}.tagged.bam" \
    isoseq refine --require-polya -j $THREADS --log-level INFO \
      "${P}.tagged.bam" "$PRIMERS" "${P}.fltnc.bam" \
    || { echo "$LIB refine" >> "$FAILED"; continue; }

  step "$LIB" correct "${P}.fltnc.bam" \
    isoseq correct --barcodes "$WL_RC" -j $THREADS --log-level INFO \
      "${P}.fltnc.bam" "${P}.corrected.bam" \
    || { echo "$LIB correct" >> "$FAILED"; continue; }

  step "$LIB" sort "${P}.corrected.bam" \
    samtools sort -t CB -@ $SORT_THREADS \
      "${P}.corrected.bam" -o "${P}.corrected.CBsorted.bam" \
    || { echo "$LIB sort" >> "$FAILED"; continue; }

  step "$LIB" dedup "${P}.corrected.CBsorted.bam" \
    isoseq groupdedup --log-level INFO \
      "${P}.corrected.CBsorted.bam" "${P}.dedup.bam" \
    || { echo "$LIB dedup" >> "$FAILED"; continue; }

  step "$LIB" align "${P}.dedup.bam" \
    pbmm2 align --preset ISOSEQ --sort -j $THREADS \
      "$REF_FA" "${P}.dedup.bam" "${P}.aligned.sorted.bam" \
    || { echo "$LIB align" >> "$FAILED"; continue; }

  echo "    Moving results to NAS"
  shopt -s nullglob
  OUTFILES=( ${P}.oriented* ${P}.tagged* ${P}.fltnc* ${P}.corrected* ${P}.dedup* ${P}.aligned* )
  shopt -u nullglob
  if rsync -ahP --remove-source-files "${OUTFILES[@]}" "${NAS_OUT}/"; then
    # remove the local copy of the input Skera file
    [ -f "${NAS_IN}/${IN}" ] && rm -f "$IN" "${IN}.pbi"
    echo "    $LIB done"
  else
    echo "    [FAILED] move to NAS, local files kept"; echo "$LIB move_to_nas" >> "$FAILED"
  fi
done

echo "Finished. Failures (if any) are in $FAILED"
