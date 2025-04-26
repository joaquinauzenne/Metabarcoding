#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Master eDNA pipeline:
#   1) Preprocess   (Cutadapt trimming → *_trimmed.fastq.gz)
#   2) QC & merge   (FastQC, NanoPlot, MultiQC, per-barcode merge)
#   3) Consensus    (NGSpeciesID @90%, + CD-HIT-EST@100% + MAFFT + consensus.pl)
#   4) Taxonomy     (blast_taxonomy.sh)
#
# Usage:
#   bash pipeline.sh <RAW_DIR> <TRIMMED_DIR> <CONS_BASE> <DB_DIR> [THREADS] [RECLUSTER]
#
#   RAW_DIR      Directory with raw .fastq(.gz) files
#   TRIMMED_DIR  Output dir for trimmed FASTQs
#   CONS_BASE    Base output dir for consensus step
#   DB_DIR       Directory containing your formatted BLAST nt DB
#   THREADS      (optional) Number of threads [default: 4]
#   RECLUSTER    (optional) true|false — whether to recluster pooled consensi
# ─────────────────────────────────────────────────────────────────────────────

RAW_DIR="${1%/}"
TRIMMED_DIR="${2%/}"
CONS_BASE="${3%/}"
DB_DIR="${4%/}"
THREADS="${5:-4}"
RECLUSTER="${6:-true}"

echo "▶ Starting pipeline:"
echo "    RAW_DIR     = $RAW_DIR"
echo "    TRIMMED_DIR = $TRIMMED_DIR"
echo "    CONS_BASE   = $CONS_BASE"
echo "    DB_DIR      = $DB_DIR"
echo "    THREADS     = $THREADS"
echo "    RECLUSTER   = $RECLUSTER"
echo

# ─── 0) Preflight Checks ─────────────────────────────────────────────────────
declare -a NEED_EXE=(
  "scripts/preprocess_reads.sh"
  "scripts/qc_nanoplot_multiqc.sh"
  "scripts/run_consensus.sh"
  "NGSpeciesID"
  "cd-hit-est"
  "mafft"
  "cutadapt"
  "blastn"
  "perl"
)
declare -a NEED_FILE=(
  "scripts/consensus.pl"
)
declare -a NEED_DIR=(
  "$RAW_DIR"
)

echo "🛫 Preflight checks..."
for exe in "${NEED_EXE[@]}"; do
  if [[ -x "$exe" ]] || command -v "$(basename "$exe")" &>/dev/null; then
    echo "✔  $exe"
  else
    echo "❌  Missing: $exe"
    exit 1
  fi
done

for file in "${NEED_FILE[@]}"; do
  if [[ -f "$file" ]]; then
    echo "✔  $file"
  else
    echo "❌  Missing file: $file"
    exit 1
  fi
done

for dir in "${NEED_DIR[@]}"; do
  if [[ -d "$dir" ]]; then
    echo "✔  DIR  $dir"
  else
    echo "❌  Missing directory: $dir"
    exit 1
  fi
done

echo "✅ Preflight OK."
echo

# ─── 1) Preprocessing ─────────────────────────────────────────────────────────
if find "$TRIMMED_DIR" -type f -name "*_trimmed.fastq.gz" | grep -q .; then
  echo "✔ SKIP Phase 1: trimmed reads already exist."
else
  echo "▶ Phase 1: preprocess reads → $TRIMMED_DIR"
  bash scripts/preprocess_reads.sh "$RAW_DIR" "$TRIMMED_DIR"
fi

# ─── 2) QC & Merge ───────────────────────────────────────────────────────────
echo
echo "▶ Phase 2: QC & merge trimmed FASTQs"
bash scripts/qc_nanoplot_multiqc.sh "$TRIMMED_DIR" "$THREADS"

# ─── 3) Clustering & Consensus ───────────────────────────────────────────────
echo
echo "▶ Phase 3: clustering & consensus (reclustering=${RECLUSTER})"
bash scripts/run_consensus.sh \
     "$TRIMMED_DIR/merged_fastq" \
     "$CONS_BASE" \
     "$THREADS" \
     "$RECLUSTER"

# ─────────────────────────────────────────────────────────────────────────────
# 4) Taxonomy
# ─────────────────────────────────────────────────────────────────────────────
echo
echo "▶ Phase 4: taxonomic assignment & reporting"
bash scripts/blast_taxonomy.sh \
     "$CONS_BASE" \
     "$DB_DIR" \
     "$THREADS"

echo
echo "✅ Pipeline complete."
echo " • Raw data:              $RAW_DIR"
echo " • Trimmed reads:         $TRIMMED_DIR"
echo " • Merged FASTQs:         $TRIMMED_DIR/merged_fastq"
echo " • NanoPlot results:      $TRIMMED_DIR/nanoplot_results_trimmed"
echo " • MultiQC report:        $TRIMMED_DIR/multiqc_report_trimmed/multiqc_report.html"
echo " • Per-barcode consensus: $CONS_BASE/per_barcode"
echo " • Reclustering:          $RECLUSTER"
echo " • Per-barcode consensus: $CONS_BASE/per_barcode/"
echo " • Final consensus FASTA: $CONS_BASE/combined/final_consensus.fasta"
echo " • Taxonomy results:      $CONS_BASE/per_barcode/*/taxonomy/, $CONS_BASE/combined/taxonomy/"
