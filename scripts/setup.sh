#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# Script : setup.sh
# Purpose: Initialize working directory of Metabarcoding and all necessary
#          scripts and data from Github repository. Ensure that NCBI nt
#          database is present and up-to-date within specified directory
#          (i.e. blastdb). Use Apptainer to activate the .sif image within
#          the Apptainer directory. Ensure that eDNA is present within the
#          Metabarcoding directory.
#
# Usage:
#   bash setup.sh <RAW_DIR> <DB_DIR> [THREADS] [SETTINGS_FILE] [MODE]
#
#   RAW_DIR       eDNA directory containing raw FASTQ(.gz) files
#   DB_DIR        Directory containing formatted BLAST nt DB (nt.*) OR a
#                 remote DB name (e.g. nt)
#   THREADS       Number of threads to use (default: all available cores)
#   SETTINGS_FILE Optional path to a settings file for job scheduling
#   MODE          Pipeline mode: otu | asv (default: otu)
#                   otu: standard NGSpeciesID clustering
#                   asv: vsearch derep + chimera removal (requires vsearch in container)
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
RAW_DIR="${1%/}"         # eDNA raw reads directory
DB_DIR="${2%/}"          # BLAST DB dir or remote name
THREADS="${3:-$(nproc)}" # number of threads (default: all cores)
SETTINGS_FILE="${4:-}"   # optional settings file for job scheduling
MODE="${5:-otu}"         # pipeline mode: otu | asv

# Validate core inputs
if [[ ! -d "$RAW_DIR" ]]; then
  printf "[ERROR] RAW_DIR '%s' does not exist\n" "$RAW_DIR" >&2
  exit 1
fi
if [[ -z "$DB_DIR" ]]; then
  printf "[ERROR] You must specify a BLAST database directory\n" >&2
  exit 1
fi

printf "[INFO] Starting setup for eDNA pipeline\n"
printf "  RAW_DIR:  %s\n" "$RAW_DIR"
printf "  DB_DIR:   %s\n" "$DB_DIR"
printf "  THREADS:  %s\n" "$THREADS"
printf "  MODE:     %s\n" "$MODE"
[[ -n "$SETTINGS_FILE" ]] && printf "  SETTINGS: %s\n" "$SETTINGS_FILE"

#─────────────────────────────────────────────────────────────────────────────
# Load environment and enter container
#─────────────────────────────────────────────────────────────────────────────
if command -v module &>/dev/null; then
  module load tacc-apptainer || true
  printf "[OK]   Loaded tacc-apptainer module\n"
  printf "[INFO] Launching Apptainer shell in edna-pipeline.sif\n"
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
  "Rscript"
)
declare -a NEED_FILE=(
  "scripts/consensus.pl"
)
declare -a NEED_DIR=(
  "$RAW_DIR"
  "$DB_DIR"
)

printf "[INFO] Preflight checks...\n"
for exe in "${NEED_EXE[@]}"; do
  if [[ -x "$exe" ]] || command -v "$(basename "$exe")" &>/dev/null; then
    printf "[OK]   %s\n" "$exe"
  else
    printf "[ERROR] Missing executable: %s\n" "$exe" >&2
    exit 1
  fi
done
for file in "${NEED_FILE[@]}"; do
  if [[ -f "$file" ]]; then
    printf "[OK]   %s\n" "$file"
  else
    printf "[ERROR] Missing file: %s\n" "$file" >&2
    exit 1
  fi
done
for dir in "${NEED_DIR[@]}"; do
  if [[ -d "$dir" ]]; then
    printf "[OK]   DIR  %s\n" "$dir"
  else
    printf "[ERROR] Missing directory: %s\n" "$dir" >&2
    exit 1
  fi
done

printf "[OK]   Preflight OK.\n"

#─────────────────────────────────────────────────────────────────────────────
# Check for required tools
#─────────────────────────────────────────────────────────────────────────────
for cmd in blastn apptainer makeblastdb update_blastdb.pl esearch efetch ktImportText sbatch; do
  command -v "$cmd" &>/dev/null || {
    printf "[ERROR] '%s' not found in PATH\n" "$cmd" >&2
    exit 1
  }
done

# ASV mode requires vsearch for dereplicate and chimera removal.
# If not present, rebuild the conda env with: conda install -c bioconda vsearch
if [[ "${MODE}" == "asv" ]]; then
  if ! command -v vsearch &>/dev/null; then
    printf "[ERROR] vsearch not found -- required for MODE=asv.\n" >&2
    printf "        Rebuild conda env: conda install -y -c bioconda vsearch\n" >&2
    printf "        Then repack and rebuild Apptainer/edna-pipeline.sif.\n" >&2
    exit 1
  fi
  printf "[OK]   vsearch available for ASV mode\n"
fi

printf "[OK]   All required tools are available\n"

#─────────────────────────────────────────────────────────────────────────────
# Ensure BLAST nt database is present or prompt download
#─────────────────────────────────────────────────────────────────────────────
if [[ -d "$DB_DIR" ]]; then
  if ! ls "$DB_DIR"/nt.* &>/dev/null; then
    read -p "BLAST nt DB not found in $DB_DIR. Download? (yes/no): " yn
    case "$yn" in
      [Yy]* )
        printf "[INFO] Downloading nt DB into %s...\n" "$DB_DIR"
        update_blastdb.pl --decompress nt --blastdb_version 5 --destination "$DB_DIR" --num_threads $THREADS
        ;;
      * )
        printf "[WARN]  Database download skipped. Exiting.\n" >&2
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

printf "[OK]   Directory structure initialized\n"

#─────────────────────────────────────────────────────────────────────────────
# Logging setup
#─────────────────────────────────────────────────────────────────────────────
LOG="logs/setup_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee "$LOG") 2>&1

printf "[INFO] Logging to %s\n" "$LOG"

#─────────────────────────────────────────────────────────────────────────────
# Error trapping
#─────────────────────────────────────────────────────────────────────────────
trap 'printf "[ERROR] Setup failed at line %s\n" "$LINENO" >&2; exit 1' ERR

printf "[OK]   Error trap installed\n"

#─────────────────────────────────────────────────────────────────────────────
# Record tool versions
#─────────────────────────────────────────────────────────────────────────────
printf "[INFO] Versions for reproducibility:\n"
printf "  blastn:      %s\n" "$(blastn -version | head -1)"
printf "  apptainer:   %s\n" "$(apptainer --version)"
printf "  makeblastdb: %s\n" "$(makeblastdb -version | head -1)"

printf "[OK]   Versions recorded\n"

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

  printf "[INFO] Job scheduler script: %s\n" "$jobScheduler"
  if [[ "$work_env" == "HPC" ]]; then
    printf "[INFO] Submitting via sbatch...\n"
    sbatch "$jobScheduler"
  else
    printf "[INFO] Running locally (sbatch not available)...\n"
    bash "$jobScheduler"
  fi
fi

printf "[OK]   Setup complete.\n"
