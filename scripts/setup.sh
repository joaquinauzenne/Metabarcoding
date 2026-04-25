#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# Script : setup.sh
# Purpose: Host-side preflight checker for the eDNA metabarcoding pipeline.
#          Run this OUTSIDE the Apptainer container, on the Lonestar6 login or
#          compute node, before submitting batch_job.sh.  It validates that:
#            - The Apptainer .sif image exists and key tools are accessible
#            - Required pipeline scripts are present
#            - The BLAST database directory and raw read directory exist
#            - Output directory structure is initialised
#
# Usage:
#   bash scripts/setup.sh <RAW_DIR> <DB_DIR> [THREADS] [SETTINGS_FILE] [MODE]
#
#   RAW_DIR       Directory containing raw FASTQ(.gz) files
#   DB_DIR        Directory containing the formatted BLAST DB (e.g. blastdb/mammalia_nt)
#   THREADS       Number of threads (default: all available cores)
#   SETTINGS_FILE Optional settings file (currently unused; reserved for future use)
#   MODE          Pipeline mode: otu | asv  (default: otu)
# ─────────────────────────────────────────────────────────────────────────────

usage() {
  grep '^#' "$0" | cut -c2-
  exit 1
}
[[ "${1:-}" =~ ^-h|--help$ ]] && usage

# ─────────────────────────────────────────────────────────────────────────────
# Arguments
# ─────────────────────────────────────────────────────────────────────────────
RAW_DIR="${1%/}"
DB_DIR="${2%/}"
THREADS="${3:-$(nproc)}"
SETTINGS_FILE="${4:-}"
MODE="${5:-otu}"

if [[ ! -d "$RAW_DIR" ]]; then
  printf "[ERROR] RAW_DIR '%s' does not exist\n" "$RAW_DIR" >&2; exit 1
fi
if [[ -z "$DB_DIR" ]]; then
  printf "[ERROR] DB_DIR must be specified\n" >&2; exit 1
fi

printf "[INFO] eDNA pipeline preflight\n"
printf "  RAW_DIR:  %s\n" "$RAW_DIR"
printf "  DB_DIR:   %s\n" "$DB_DIR"
printf "  THREADS:  %s\n" "$THREADS"
printf "  MODE:     %s\n" "$MODE"

# ─────────────────────────────────────────────────────────────────────────────
# 1. Apptainer — load module if available, then verify .sif exists
# ─────────────────────────────────────────────────────────────────────────────
SIF="Apptainer/edna-pipeline.sif"

# Note: load tacc-apptainer before running this script:
#   module load tacc-apptainer
if ! command -v apptainer &>/dev/null; then
  printf "[ERROR] apptainer not found in PATH. Run: module load tacc-apptainer\n" >&2
  exit 1
fi
printf "[OK]   apptainer: %s\n" "$(apptainer --version)"

if [[ ! -f "$SIF" ]]; then
  printf "[ERROR] Container image not found: %s\n" "$SIF" >&2
  exit 1
fi
printf "[OK]   Container image: %s\n" "$SIF"

# ─────────────────────────────────────────────────────────────────────────────
# 2. Pipeline scripts
# ─────────────────────────────────────────────────────────────────────────────
printf "[INFO] Checking pipeline scripts...\n"
for f in scripts/pipeline.sh scripts/setup.sh scripts/batch_job.sh \
          scripts/consensus.pl scripts/lca_assignment.R; do
  if [[ -f "$f" ]]; then
    printf "[OK]   %s\n" "$f"
  else
    printf "[ERROR] Missing: %s\n" "$f" >&2; exit 1
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 3. Tool check inside container
# ─────────────────────────────────────────────────────────────────────────────
printf "[INFO] Checking tools inside container...\n"
TOOLS=(cutadapt fastqc NanoPlot multiqc blastn NGSpeciesID mafft
       cd-hit-est seqtk perl Rscript ktImportText)
[[ "$MODE" == "asv" ]] && TOOLS+=(vsearch)

for t in "${TOOLS[@]}"; do
  if apptainer exec "$SIF" bash -c "command -v '$t'" &>/dev/null; then
    printf "[OK]   %s\n" "$t"
  else
    printf "[ERROR] Tool not found in container: %s\n" "$t" >&2; exit 1
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 4. BLAST database
# DB_DIR is a BLAST path prefix (e.g. blastdb/nt_euk), not a directory.
# blastn is called as: blastn -db "$DB_DIR"
# ─────────────────────────────────────────────────────────────────────────────
printf "[INFO] Checking BLAST database prefix: %s\n" "$DB_DIR"
DB_PARENT=$(dirname "$DB_DIR")
if [[ ! -d "$DB_PARENT" ]]; then
  printf "[ERROR] DB parent directory '%s' does not exist\n" "$DB_PARENT" >&2; exit 1
fi
if ls "${DB_DIR}".*.nhr &>/dev/null || ls "${DB_DIR}".nhr &>/dev/null; then
  printf "[OK]   BLAST DB index files found for prefix: %s\n" "$DB_DIR"
else
  printf "[WARN] No .nhr index files found matching '%s.*'\n" "$DB_DIR" >&2
  printf "       Available prefixes in %s:\n" "$DB_PARENT" >&2
  ls "$DB_PARENT"/*.nal 2>/dev/null | sed 's/\.nal$//' | xargs -I{} basename {} | sort -u | sed 's/^/         /' >&2 || true
  printf "       Run: DOWNLOAD_DB=true sbatch scripts/mammalia_nt_db.sh\n" >&2
fi

# ─────────────────────────────────────────────────────────────────────────────
# 5. Raw reads
# ─────────────────────────────────────────────────────────────────────────────
printf "[INFO] Checking raw reads in %s...\n" "$RAW_DIR"
if ls "$RAW_DIR"/**/*.fastq* &>/dev/null || ls "$RAW_DIR"/*.fastq* &>/dev/null; then
  printf "[OK]   FASTQ files found\n"
else
  printf "[WARN] No FASTQ files found in %s (may be in barcode subdirectories)\n" "$RAW_DIR"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 6. Output directory structure
# ─────────────────────────────────────────────────────────────────────────────
printf "[INFO] Initialising output directories...\n"
mkdir -p logs
printf "[OK]   logs/ ready\n"

# ─────────────────────────────────────────────────────────────────────────────
# 7. batch_job.sh allocation check
# ─────────────────────────────────────────────────────────────────────────────
if grep -q '<your-allocation>' scripts/batch_job.sh 2>/dev/null; then
  printf "[WARN] batch_job.sh still has placeholder '<your-allocation>' — fill in before sbatch\n"
else
  ALLOC=$(grep '#SBATCH -A' scripts/batch_job.sh | awk '{print $3}')
  printf "[OK]   batch_job.sh allocation: %s\n" "$ALLOC"
fi

printf "[OK]   Setup complete.\n"
