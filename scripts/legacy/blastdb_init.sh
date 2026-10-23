#!/usr/bin/env bash
#
# DEPRECATED — 2026-04-25
# ---------------------------------------------------------------
# This script downloads the full NCBI `nt` database (~600 GB) via
# update_blastdb.pl. The active DB-build path is now:
#
#     DOWNLOAD_DB=true sbatch scripts/mammalia_nt_db.sh
#
# which fetches the smaller `nt_euk` (Eukaryota only, ~50 GB) via
# rsync and extracts a Mammalia-only subset for the eDNA pipeline.
# update_blastdb.pl was not chosen for the active path because it
# hard-codes Perl Net::FTP (port 21) and TACC's firewall blocks
# FTP — see contextLiterature/eln-blastdb-init-deprecation-2026-04-25.md.
#
# Kept here for historical reference and as a starting point if
# full-nt is ever needed again. Not validated against the post-
# 2026-04 module/XALT fixes — re-test before using.
# ---------------------------------------------------------------
#
#SBATCH --job-name=blastdb_init
#SBATCH --output=blastdb_init.o%j
#SBATCH --error=blastdb_init.e%j
#SBATCH --partition=development
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00

####— Fail-fast and error trap —####
set -o errexit    # Exit on non-zero :contentReference[oaicite:7]{index=7}
set -o pipefail   # Fail on any pipe component :contentReference[oaicite:9]{index=9}
set -o errtrace   # ERR trap in subshells :contentReference[oaicite:10]{index=10}
# nounset is enabled AFTER module load (see below) — tacc-apptainer's
# bash-completion references BASH_COMPLETION_DEBUG without a default,
# which aborts the job under `set -u` during module load.

trap 'echo "[ERROR] Line $LINENO: \`$BASH_COMMAND\` failed with status $?" >&2' ERR

####— 1) Modules & working directory —####
echo "[INFO] Loading modules and changing to Apptainer dir"
module load tacc-apptainer                         # Load Apptainer :contentReference[oaicite:12]{index=12}

set -o nounset    # Error on undefined vars :contentReference[oaicite:8]{index=8}

cd "$WORK/Metabarcoding/Apptainer" || { 
  echo "[ERROR] Cannot cd to Apptainer dir $WORK/Metabarcoding/Apptainer" >&2
  exit 1
}

####— 2) Define & stage directories —####
echo "[INFO] Setting up BLASTDB dirs"
export BLASTDB="$WORK/Metabarcoding/blastdb/nt"
mkdir -p "$BLASTDB" || { echo "[ERROR] mkdir $BLASTDB failed"; exit 1; }

export SCRATCH_DB="$SCRATCH/Metabarcoding/blastdb/nt"
mkdir -p "$SCRATCH_DB" || { echo "[ERROR] mkdir $SCRATCH_DB failed"; exit 1; }

echo "[INFO] Syncing existing DB from $BLASTDB to $SCRATCH_DB"
rsync -a "$BLASTDB/" "$SCRATCH_DB/" || { 
  echo "[ERROR] rsync BLASTDB -> SCRATCH failed"; exit 1; 
}
cd "$SCRATCH_DB" || { echo "[ERROR] cd to $SCRATCH_DB failed"; exit 1; }

####— 3) Download & decompress nt —####
echo "[INFO] Verifying update_blastdb.pl is available"
command -v update_blastdb.pl >/dev/null || {
  echo "[ERROR] update_blastdb.pl not found in PATH"; exit 1;
}

echo "[INFO] Downloading NCBI nt database"
if ! apptainer exec \
     --bind "$SCRATCH":"$SCRATCH" \
     --bind "$WORK":"$WORK" \
     edna-pipeline.sif \
     bash -lc 'update_blastdb.pl --decompress nt --num_threads 16'; then
  echo "[ERROR] Failed to decompress nt database"; exit 1;
fi

####— 4) (Optional) taxonomy data —####
echo "[INFO] Downloading taxonomy database"
if ! apptainer exec \
     --bind "$SCRATCH":"$SCRATCH" \
     --bind "$WORK":"$WORK" \
     edna-pipeline.sif \
     bash -lc 'update_blastdb.pl --decompress taxdb'; then
  echo "[ERROR] Failed to decompress taxdb"; exit 1;
fi

gunzip taxdb.tar.gz || { echo "[ERROR] gunzip taxdb.tar.gz failed"; exit 1; }

####— 5) Sync back to WORK —####
echo "[INFO] Syncing back to $BLASTDB"
rsync -a "$SCRATCH_DB/" "$BLASTDB/" || { 
  echo "[ERROR] rsync SCRATCH -> WORK failed"; exit 1; 
}

echo "[SUCCESS] BLAST nt database initialized in $BLASTDB"
