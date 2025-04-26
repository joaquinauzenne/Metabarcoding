#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Script : preprocess_reads.sh
# Purpose: Phase 1 of the eDNA pipeline — (optional) adapter trimming +
#          poly‑A tails + quality- and length-filtering using Cutadapt.
#
# Usage:
#    bash preprocess_reads.sh /path/to/input_fastq_dir /path/to/output_trimmed_dir
# ─────────────────────────────────────────────────────────────────────────────

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# User‑tweakable parameters
MIN_QUAL=12            # minimum Phred quality cutoff for end‑trimming
MIN_LEN=10             # minimum read length to keep
MAX_LEN=1100           # maximum read length to keep

mkdir -p "$OUTPUT_DIR"

echo "▶ Phase 1: preprocessing reads"
echo "   Quality cutoff: $MIN_QUAL"
echo "   Length range : $MIN_LEN–$MAX_LEN bp"

# Find both .fastq and .fastq.gz
find "$INPUT_DIR" -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) | while IFS= read -r IN; do
  # strip INPUT_DIR/, then strip .fastq or .fastq.gz
  REL="${IN#"$INPUT_DIR"/}"
  base="${REL%.fastq.gz}"
  base="${base%.fastq}"

  OUT="$OUTPUT_DIR/${base}_trimmed.fastq.gz"
  mkdir -p "$(dirname "$OUT")"

  # Skip already‑processed files
  if [[ -s "$OUT" ]]; then
    echo "✔ SKIP: $REL (already exists)"
    continue
  fi

  echo "⟳ Trimming & filtering: $REL"

  cutadapt \
    --poly-a \
    --trim-n \
    --max-n 0.05 \
    --minimum-length "$MIN_LEN" \
    --maximum-length "$MAX_LEN" \
    --quality-cutoff "$MIN_QUAL" \
    --cores 0 \
    -o "$OUT" \
    "$IN"

  if [[ -s "$OUT" ]]; then
    echo "✔ Trimmed: ${base}_trimmed.fastq.gz"
  else
    echo "✖ ERROR: no output for $REL" >&2
    rm -f "$OUT"
  fi
done

echo "✔ Phase 1 complete."
