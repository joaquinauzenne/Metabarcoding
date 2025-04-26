#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# Script: sequential_cutadapt.sh
# Purpose: Sequentially trim sequence adapters from Oxford Nanopore Technologies (ONT) kits and poly‑A/T tails from all FASTQ
#          files under a directory tree using Cutadapt.
#
# Usage:
#   bash sequential_cutadapt.sh /path/to/_unzipped /path/to/output_trimmed
# ----------------------------------------------------------------------

INPUT_DIR="$"
OUTPUT_DIR="$2"

# Temp subdirectories
TMP1="$OUTPUT_DIR/step1_remove_FWD"
TMP2="$OUTPUT_DIR/step2_remove_REV"

# Create top‑level dirs
mkdir -p "$INPUT_DIR" "$TMP1" "$TMP2" "$OUTPUT_DIR"

# Loop through every .fastq
find "$INPUT_DIR" -type f -name '*.fastq' | while read -r FILE; do
  # Strip off the INPUT_DIR/ prefix to get a truly relative path
  REL_PATH="${FILE#${INPUT_DIR}/}"    # e.g. "Qia25/Barcode08/… .fastq"
  REL_DIR="$(dirname "$REL_PATH")"    # e.g. "Qia25/Barcode08"
  BASE="$(basename "$REL_PATH" .fastq)" 

  # Ensure subdirs exist
  mkdir -p "$TMP1/$REL_DIR" "$TMP2/$REL_DIR" "$OUTPUT_DIR/$REL_DIR"

  # 1) Trim forward primer (5' end)
  cutadapt \
    -g ^"$FWD_PRIMER" \
    -o "$TMP1/$REL_DIR/${BASE}_noFWD.fastq" \
    "$FILE"

  # 2) Trim reverse primer (3' end)
  cutadapt \
    -a "$REV_PRIMER"\$ \
    -o "$TMP2/$REL_DIR/${BASE}_noREV.fastq" \
    "$TMP1/$REL_DIR/${BASE}_noFWD.fastq"

  # 3) Trim homopolymer tails (>=10 A's at 3' )
  cutadapt \
    --poly-a \
    --minimum-length 60 \
    --quality-cutoff 10 \
    -o "$OUTPUT_DIR/$REL_DIR/${BASE}_trimmed.fastq" \
    "$TMP2/$REL_DIR/${BASE}_noREV.fastq"

  echo "✔ Trimmed: $REL_PATH → $REL_DIR/${BASE}_trimmed.fastq"

    # --- CLEANUP: remove intermediate files ---
    rm -f "$INTER1" "$INTER2"

done
