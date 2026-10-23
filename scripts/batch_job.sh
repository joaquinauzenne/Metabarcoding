#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Script : batch_job.sh
# Purpose: SLURM submission script for a full 12-barcode OTU run on Lonestar6.
#
# Usage (run directly — prompts for allocation + phase override):
#   bash scripts/batch_job.sh
#
# Resume from a specific phase:
#   bash scripts/batch_job.sh     # then enter e.g. "Phase3,Phase4" at prompt
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH -J edna-otu
#SBATCH -o logs/edna_%j.out
#SBATCH -e logs/edna_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH -t 08:00:00
# Note: --account is NOT hardcoded here. The script prompts for it
#       interactively when run directly, then self-submits via sbatch.

# ─── Self-submission block (runs only outside of a SLURM job) ────────────────
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf "╔══════════════════════════════════════════════════════════════╗\n"
  printf "║  batch_job.sh — eDNA pipeline SLURM submission              ║\n"
  printf "╚══════════════════════════════════════════════════════════════╝\n"

  # Show available allocations to help the user choose
  if command -v show_allocations &>/dev/null; then
    printf "\nYour available allocations:\n"
    show_allocations 2>/dev/null || true
    printf "\n"
  fi

  # Prompt 1: allocation account
  read -rp "Enter your TACC allocation account: " ACCOUNT
  if [[ -z "$ACCOUNT" ]]; then
    printf "[ERROR] No allocation entered. Aborting.\n" >&2
    exit 1
  fi

  # Prompt 2: phase selection
  printf "\nPipeline phases to run:\n"
  printf "  all              — run all phases (default)\n"
  printf "  Phase1           — preprocess only\n"
  printf "  Phase3,Phase4    — consensus + taxonomy (resume after QC)\n"
  read -rp "Phases [default=all]: " PHASE_INPUT
  MODULES="${PHASE_INPUT:-all}"
  printf "[INFO] MODULES=%s\n" "$MODULES"

  mkdir -p logs

  JOB_ID=$(MODULES="$MODULES" \
    sbatch --account="$ACCOUNT" "$(realpath "$0")" | awk '{print $NF}')

  printf "[OK]   Submitted job %s with account %s\n" "$JOB_ID" "$ACCOUNT"
  printf "       Monitor : squeue -j %s\n" "$JOB_ID"
  printf "       Log     : logs/edna_%s.out\n" "$JOB_ID"
  exit 0
fi

# ─── Everything below runs inside the SLURM job ──────────────────────────────

set -eo pipefail

# TACC's tacc-apptainer bash-completion script references BASH_COMPLETION_DEBUG
# without a default; under `set -u` this aborts the job during module load.
module load tacc-apptainer

set -u

# MODULES is passed via the self-submission block above, or can be overridden:
#   MODULES="Phase3,Phase4" bash scripts/batch_job.sh
MODULES="${MODULES:-all}"

apptainer exec Apptainer/edna-pipeline.sif \
  bash scripts/pipeline.sh \
  eDNA/2024 \
  eDNA/otu_trimmed \
  eDNA/otu \
  blastdb \
  8 true misc/primers.txt otu "$MODULES"
# DB_DIR is the dir holding the BLAST DB; the DB name comes from config
# (blast.db_name, default mammalia_nt) -> pipeline queries blastdb/mammalia_nt.
# Build it with scripts/mammalia_nt_db.sh.
