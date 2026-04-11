#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Script : batch_job.sh
# Purpose: SLURM submission script for a full 12-barcode OTU run on Lonestar6.
#          Edit <your-allocation> before submitting.
#
# Usage:
#   sbatch scripts/batch_job.sh
#
# Resume from a specific phase (e.g. Phase3 onwards):
#   sbatch --export=ALL,MODULES="Phase3,Phase4" scripts/batch_job.sh
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH -J edna-otu
#SBATCH -o logs/edna_%j.out
#SBATCH -e logs/edna_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH -t 08:00:00
#SBATCH -A <your-allocation>

set -euo pipefail

module load tacc-apptainer

# MODULES can be overridden at sbatch time via --export=ALL,MODULES="Phase3,Phase4"
MODULES="${MODULES:-all}"

apptainer exec Apptainer/edna-pipeline.sif \
  bash scripts/pipeline.sh \
  eDNA/2024 \
  eDNA/otu_trimmed \
  eDNA/otu \
  blastdb/mammalia_nt \
  8 true misc/primers.txt otu "$MODULES"
