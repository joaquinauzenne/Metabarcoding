#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Script: qc_nanoplot_multiqc.sh
# Purpose: Phase 2 QC & merge for the eDNA pipeline:
#   • FastQC on trimmed FASTQs (skip if already present)
#   • Merge *_trimmed.fastq(.gz) by Barcode → all merged FASTQs in one dir
#   • NanoPlot on merged FASTQs
#   • MultiQC report aggregation
#
# Usage:
#   bash qc_nanoplot_multiqc.sh <TRIMMED_DIR> [THREADS]
#
# Example:
#   bash qc_nanoplot_multiqc.sh /data/eDNA/trimmed 8
# ─────────────────────────────────────────────────────────────────────────────

TRIMMED_DIR="${1:?Usage: $0 TRIMMED_DIR [THREADS]}"
THREADS="${2:-4}"

FASTQC_OUT="$TRIMMED_DIR/fastqc_results_trimmed"
MERGED_DIR="$TRIMMED_DIR/merged_fastq"
NANOPLOT_OUT="$TRIMMED_DIR/nanoplot_results_trimmed"
MULTIQC_OUT="$TRIMMED_DIR/multiqc_report_trimmed"

# ─── 0) Prep output dirs ────────────────────────────────────────────────────
mkdir -p \
  "$FASTQC_OUT" \
  "$MERGED_DIR" \
  "$NANOPLOT_OUT" \
  "$MULTIQC_OUT"

# ─── 1) FastQC (skip if already done) ────────────────────────────────────────
if [[ -d "$FASTQC_OUT" && -n "$(ls -A "$FASTQC_OUT")" ]]; then
  echo "✅ FastQC outputs already present; skipping FastQC."
else
  echo "🔍 Running FastQC on all *_trimmed.fastq(.gz)..."
  find "$TRIMMED_DIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) | \
  while IFS= read -r fq; do
    echo "▶ DEBUG: Starting FastQC on $fq"
    if fastqc --threads "$THREADS" --outdir "$FASTQC_OUT" "$fq"; then
      echo "✔ DEBUG: FastQC succeeded for $fq"
    else
      echo "✖ DEBUG: FastQC FAILED for $fq" >&2
    fi
  done
fi

# ─── 2) Merge trimmed FASTQs by Barcode into one flat directory ─────────────
echo "🔗 Merging trimmed FASTQs by Barcode into $MERGED_DIR..."
# Find each directory containing any *_trimmed.fastq(.gz), then merge its contents
find "$TRIMMED_DIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) \
  -printf '%h\n' | sort -u | while IFS= read -r BDIR; do
  # create a safe sample tag, replacing "/" with "_"
  SAMPLE_TAG="${BDIR#${TRIMMED_DIR}/}"
  SAMPLE_TAG="${SAMPLE_TAG//\//_}"
  MERGED_FASTQ="$MERGED_DIR/${SAMPLE_TAG}_merged.fastq"

  echo "  • Merging $BDIR/*_trimmed.fastq* → $MERGED_FASTQ"
  {
    find "$BDIR" -type f \( -name "*_trimmed.fastq" -o -name "*_trimmed.fastq.gz" \) \
      | sort | \
      while IFS= read -r in; do
        if [[ "$in" == *.gz ]]; then
          zcat "$in"
        else
          cat "$in"
        fi
      done
  } > "$MERGED_FASTQ"
done

# ─── 3) NanoPlot on merged FASTQs ──────────────────────────────────────────
echo "📊 Running NanoPlot on merged FASTQs..."
find "$MERGED_DIR" -type f -name "*_merged.fastq" | \
while IFS= read -r mfq; do
  SAMPLE="$(basename "$mfq" _merged.fastq)"
  echo "▶ DEBUG: NanoPlot on $mfq"
  if NanoPlot \
      --fastq "$mfq" \
      --outdir "$NANOPLOT_OUT/$SAMPLE" \
      --threads "$THREADS" \
      --N50 \
      --plots hex dot kde; then
    echo "✔ DEBUG: NanoPlot complete for $mfq"
  else
    echo "✖ DEBUG: NanoPlot FAILED for $mfq" >&2
  fi
done

# ─── 4) MultiQC aggregation of FastQC ──────────────────────────────────────
echo "📈 Aggregating FastQC reports with MultiQC..."
if multiqc "$FASTQC_OUT" --outdir "$MULTIQC_OUT"; then
  echo "✔ MultiQC report generated at $MULTIQC_OUT/multiqc_report.html"
else
  echo "✖ MultiQC failed" >&2
fi

echo "✅ QC & merge complete."
echo "   FastQC:    $FASTQC_OUT"
echo "   Merged:    $MERGED_DIR"
echo "   NanoPlot:  $NANOPLOT_OUT"
echo "   MultiQC:   $MULTIQC_OUT/multiqc_report.html"
