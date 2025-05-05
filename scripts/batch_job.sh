#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Script : batch_job.sh
# Purpose: Makes a batch job in TACC to run a full soil sample analysis
#
# Usage:
#   Log into TACC, then run this script to process pipeline.sh
#   bash setup.sh <CONS_BASE> <DB_DIR> [threads]
#
#   CONS_BASE:  base dir of Phase 3 output (contains per_barcode/ and combined/)
#   DB_DIR:     directory containing formatted BLAST nt DB (nt.*) OR a
#               remote DB name (e.g. nt, nr, refseq_rna)
#   threads:    number of threads to use (default: 4)
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH -J myjob           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 1               # Total # of mpi tasks
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH -A myproject       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

