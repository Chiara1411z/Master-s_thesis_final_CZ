#!/bin/bash
set -uo pipefail

# Input and Output
# Input are HiFi files
HIFI_NAS="/.../"
ADAPTERS_NAS=".../"
# OUTPUT (NAS)
SKERA_NAS="/.../" 
QC_NAS="/.../"  
# Local working directory
WORKDIR="/.../"
# -----------------------------------------------------------------------------
THREADS=8
SKERA_THREADS=16

declare -A MOVIE=(
...
)
declare -A SKERA=(
...
)

if [ $# -gt 0 ]; then LIBS=("$@"); else LIBS=(...); fi

mkdir -p "$WORKDIR" "$SKERA_NAS" "$QC_NAS"
cd "$WORKDIR" || { echo "Cannot enter $WORKDIR"; exit 1; }
FAILED="${QC_NAS}/failed_qc.log"; : > "$FAILED"
YIELD="${QC_NAS}/yield_summary.tsv"
[ -f "$YIELD" ] || printf "library\thifi_reads\ts_reads\tyield_multiplier\tmean_s_read_len\n" > "$YIELD"

rsync -ah "$ADAPTERS_NAS" . || { echo "Adapter fasta not found: $ADAPTERS_NAS"; exit 1; }
ADAPTERS="$(basename "$ADAPTERS_NAS")"

# move_to_nas 
move_to_nas() {
  local dest=$1; shift
  mkdir -p "$dest"
  rsync -ahP --remove-source-files --ignore-existing --ignore-missing-args "$@" "$dest/" || return 1
  for item in "$@"; do [ -d "$item" ] && find "$item" -type d -empty -delete; done
  return 0
}

for LIB in "${LIBS[@]}"; do
  M="${MOVIE[$LIB]:-}"
  if [ -z "$M" ]; then echo "Unknown library $LIB, skipping"; continue; fi
  HIFI="${M}.hifi_reads.bam"
  SK="${SKERA[$LIB]}"

  if grep -qP "^${LIB}\t" "$YIELD"; then echo "$LIB already in $YIELD, skipping"; continue; fi

  echo "=========================================="
  echo " $LIB  ($M)"
  echo "=========================================="

  # 1. NAS -> local: HiFi bam + pbi
  echo "--> Copying HiFi files from NAS"
  rsync -ahP "${HIFI_NAS}/${M}.hifi_reads."* . || { echo "$LIB copy_hifi" >> "$FAILED"; continue; }

  # 2. NanoPlot on raw HiFi
  echo "--> NanoPlot HiFi"
  NanoPlot --ubam "$HIFI" --threads $THREADS -o "nanoplot_hifi_${LIB}/" \
    || { echo "$LIB nanoplot_hifi" >> "$FAILED"; continue; }
  cat "nanoplot_hifi_${LIB}/NanoStats.txt"

  # 3. Skera: reuse the existing file from the NAS, otherwise create a new one
  NEW_SKERA=0
  if [ -f "${SKERA_NAS}/${SK}" ]; then
    echo "--> ${SK} already exists on NAS, copying it in (not recreating)"
    rsync -ahP --ignore-missing-args "${SKERA_NAS}/${SK}" "${SKERA_NAS}/${SK}.pbi" . \
      || { echo "$LIB copy_skera" >> "$FAILED"; continue; }
  else
    echo "--> skera split -> $SK"
    skera split -j $SKERA_THREADS --log-level INFO "$HIFI" "$ADAPTERS" "$SK" \
      > "skera_split_${LIB}.log" 2>&1 \
      || { echo "$LIB skera_split (HiFi file incomplete?)" >> "$FAILED"; rm -f "$SK" "${SK}.pbi"; continue; }
    NEW_SKERA=1
  fi

  # 4. S-read QC
  echo "--> NanoPlot S-reads + samtools stats"
  NanoPlot --ubam "$SK" --threads $THREADS -o "nanoplot_sreads_${LIB}/" \
    || { echo "$LIB nanoplot_sreads" >> "$FAILED"; continue; }
  samtools stats -@ $THREADS "$SK" > "skera_${LIB}_stats.txt" \
    || { echo "$LIB samtools_stats" >> "$FAILED"; continue; }
  grep "^SN" "skera_${LIB}_stats.txt" | grep -E "sequences:|average length:|total length:"

  # 5. Yield multiplier = S-reads / parent HiFi reads
  HIFI_N=$(samtools view -c -@ $THREADS "$HIFI")
  S_N=$(awk -F'\t' '$1=="SN" && $2=="sequences:" {print $3}' "skera_${LIB}_stats.txt")
  S_LEN=$(awk -F'\t' '$1=="SN" && $2=="average length:" {print $3}' "skera_${LIB}_stats.txt")
  MULT=$(awk -v s="$S_N" -v h="$HIFI_N" 'BEGIN{printf "%.2f", s/h}')
  echo "--> HiFi reads: $HIFI_N | S-reads: $S_N | multiplier: $MULT | mean S-read length: $S_LEN"
  printf "%s\t%s\t%s\t%s\t%s\n" "$LIB" "$HIFI_N" "$S_N" "$MULT" "$S_LEN" >> "$YIELD"

  # 6. Local -> NAS, then clean up local copies
  echo "--> Moving results to NAS"
  move_to_nas "$QC_NAS" "nanoplot_hifi_${LIB}" "nanoplot_sreads_${LIB}" "skera_${LIB}_stats.txt" "skera_split_${LIB}.log" \
    || echo "$LIB move_qc" >> "$FAILED"
  if [ $NEW_SKERA -eq 1 ]; then
    move_to_nas "$SKERA_NAS" "$SK" "${SK}.pbi" || echo "$LIB move_skera" >> "$FAILED"
  else
    rm -f "$SK" "${SK}.pbi" 
  fi
  rm -f "$HIFI" "${HIFI}.pbi" "${M}.hifi_reads."*
  echo "    $LIB done"
done

# 7. Aggregate everything on the NAS
multiqc "$QC_NAS" -o "${QC_NAS}/multiqc_report" -f

echo "Done. Yield table: $YIELD | failures (if any): $FAILED"
