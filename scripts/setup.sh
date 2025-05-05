#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# Script : setup.sh
# Purpose: Initialize working directory of Metabarcoding and all necessary scripts and data from Github repository.
#          Ensure that NCBI nt database is present and up-to-date within specified directory (i.e. blastdb).
#          Use Apptainer to activate the .sif image within the Apptainer directory.
#          Ensure that eDNA is present within the Metabarcoding directory.
#
# Usage:
#   bash setup.sh <RAW_DIR> <DB_DIR> [threads]
#
#   RAW_DIR:   eDNA directory containing raw FASTQ(.gz) files
#   DB_DIR:    directory containing formatted BLAST nt DB (nt.*) OR a
#              remote DB name (e.g. nt)
#   threads:   number of threads to use (default: all available cores)
# ─────────────────────────────────────────────────────────────────────────────

#─────────────────────────────────────────────────────────────────────────────
# Print help message
#─────────────────────────────────────────────────────────────────────────────
usage() {
  grep '^#' "$0" | cut -c2-
  exit 1
}
[[ "${1:-}" =~ ^-h|--help$ ]] && usage

#─────────────────────────────────────────────────────────────────────────────
# Parse and validate arguments
#─────────────────────────────────────────────────────────────────────────────
RAW_DIR="${1%/}"        # eDNA raw reads directory
DB_DIR="${2%/}"         # BLAST DB dir or remote name
THREADS="${3:-$(nproc)}" # number of threads (default: all cores)
SETTINGS_FILE="${4:-}"   # optional settings file for job scheduling

# Validate core inputs
if [[ ! -d "$RAW_DIR" ]]; then
  echo "ERROR: RAW_DIR directory '$RAW_DIR' does not exist"
  exit 1
fi
if [[ -z "$DB_DIR" ]]; then
  echo "ERROR: You must specify a BLAST database directory"
  exit 1
fi

echo "▶ Starting setup for eDNA pipeline"
echo "  • Raw directory  : $RAW_DIR"
echo "  • DB directory   : $DB_DIR"
echo "  • Threads        : $THREADS"
[[ -n "$SETTINGS_FILE" ]] && echo "  • Settings file  : $SETTINGS_FILE"

#─────────────────────────────────────────────────────────────────────────────
# Load environment and enter container
#─────────────────────────────────────────────────────────────────────────────
if command -v module &>/dev/null; then
  module load tacc-apptainer || true
  echo "✔ Loaded tacc-apptainer module"
  echo "▶ Launching Apptainer shell in edna-pipeline.sif"
  apptainer shell edna-pipeline.sif
else
  source activate edna-pipeline || true
fi

#─────────────────────────────────────────────────────────────────────────────
# Preflight Checks: executables, files, directories
#─────────────────────────────────────────────────────────────────────────────
declare -a NEED_EXE=(
  "scripts/pipeline.sh"
  "scripts/preprocess_reads.sh"
  "scripts/qc_nanoplot_multiqc.sh"
  "scripts/run_consensus.sh"
  "scripts/blast_taxonomy.sh"
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
  "$DB_DIR"
)

echo "🛫 Preflight checks..."
for exe in "${NEED_EXE[@]}"; do
  if [[ -x "$exe" ]] || command -v "$(basename "$exe")" &>/dev/null; then
    echo "✔  $exe"
  else
    echo "❌  Missing executable: $exe"
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

#─────────────────────────────────────────────────────────────────────────────
# Check for required tools
#─────────────────────────────────────────────────────────────────────────────
for cmd in blastn apptainer makeblastdb update_blastdb.pl esearch efetch ktImportText sbatch; do
  command -v "$cmd" &>/dev/null || {
    echo "ERROR: '$cmd' not found in PATH" >&2
    exit 1
  }
done

echo "✔ All required tools are available"

#─────────────────────────────────────────────────────────────────────────────
# Ensure BLAST nt database is present or prompt download
#─────────────────────────────────────────────────────────────────────────────
if [[ -d "$DB_DIR" ]]; then
  if ! ls "$DB_DIR"/nt.* &>/dev/null; then
    read -p "BLAST nt DB not found in $DB_DIR. Download? (yes/no): " yn
    case "$yn" in
      [Yy]* )
        echo "Downloading nt DB into $DB_DIR..."
        update_blastdb.pl --decompress nt --blastdb_version 5 --destination "$DB_DIR" --num_threads $THREADS
        ;;
      * )
        echo "Database download skipped. Exiting." >&2
        exit 1
        ;;
    esac
  fi
fi

#─────────────────────────────────────────────────────────────────────────────
# Create output directory structure
#─────────────────────────────────────────────────────────────────────────────
mkdir -p \
  "${RAW_DIR}/per_barcode" \
  "${RAW_DIR}/combined/taxonomy" \
  logs

echo "✔ Directory structure initialized"

#─────────────────────────────────────────────────────────────────────────────
# Logging setup
#─────────────────────────────────────────────────────────────────────────────
LOG="logs/setup_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee "$LOG") 2>&1

echo "Logging to $LOG"

#─────────────────────────────────────────────────────────────────────────────
# Error trapping
#─────────────────────────────────────────────────────────────────────────────
trap 'echo "ERROR: Setup failed at line $LINENO" >&2; exit 1' ERR

echo "✔ Error trap installed"

#─────────────────────────────────────────────────────────────────────────────
# Record tool versions
#─────────────────────────────────────────────────────────────────────────────
echo "Versions for reproducibility:"
echo "  blastn:      $(blastn -version | head -1)"
echo "  apptainer:   $(apptainer --version)"
echo "  makeblastdb: $(makeblastdb -version | head -1)"

echo "✔ Versions recorded"

#─────────────────────────────────────────────────────────────────────────────
# Job scheduling function (if settings file provided)
#─────────────────────────────────────────────────────────────────────────────
if [[ -n "$SETTINGS_FILE" && -f "$SETTINGS_FILE" ]]; then
  # Parse settings
  dos2unix "$SETTINGS_FILE"
  RunAnalysisScript="scripts/pipeline.sh"
  RunName=$(grep '^RunName=' "$SETTINGS_FILE" | cut -d'=' -f2 | tr -d '"')
  OutDir=$(grep '^OutDir=' "$SETTINGS_FILE" | cut -d'=' -f2 | tr -d '"')
  GPU=$(grep '^GPU=' "$SETTINGS_FILE" | cut -d'=' -f2 | tr -d '"')
  MAIL_USER=$(grep '^mail_user=' "$SETTINGS_FILE" | cut -d'=' -f2 | tr -d '"')
  mkdir -p "$OutDir/$RunName/job_scheduler/LOGS"
  work_env='HPC'
  if ! command -v sinfo &>/dev/null; then work_env='local'; fi

  jobScheduler="$OutDir/$RunName/job_scheduler/jobScheduler.${RunName}.sh"
  cat > "$jobScheduler" <<EOF
#!/bin/bash
#SBATCH --comment=edna_pipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem-per-cpu=8000
#SBATCH --time=200:00:00
#SBATCH --output=$OutDir/$RunName/job_scheduler/LOGS/output_%j.txt
#SBATCH --error=$OutDir/$RunName/job_scheduler/LOGS/error_%j.txt
#SBATCH --job-name=${RunName}
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$MAIL_USER

EOF
