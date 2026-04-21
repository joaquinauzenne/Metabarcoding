#!/usr/bin/env bash
#SBATCH --job-name=mammalia_nt_db
#SBATCH --output=logs/mammalia_nt_%j.out
#SBATCH --error=logs/mammalia_nt_%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
# Note: --account is NOT hardcoded here. The script prompts for it
#       interactively when run directly, then self-submits via sbatch.

# ─────────────────────────────────────────────────────────────────────────────
# mammalia_nt_db.sh
# Build a Mammalia-only BLAST nucleotide database from nt_euk.
#
# Download policy:  NEVER downloads unless the user explicitly opts in.
#   - nt_euk present  → uses it as-is (warns if stale)
#   - nt_euk missing  → ERROR with instructions to opt in
#   - DOWNLOAD_DB=true → prompted interactively, or set before calling
#
# Usage (run directly — prompts for allocation + download choice):
#   bash scripts/mammalia_nt_db.sh
#
# Output:
#   $WORK/Metabarcoding/blastdb/mammalia_nt.*  (BLAST DB prefix)
#
# Taxids:
#   Mammalia : 40674  (all mammals)
#   Hominidae: 9604   (excluded — reduces human contamination signal)
# ─────────────────────────────────────────────────────────────────────────────

# ─── Self-submission block (runs only outside of a SLURM job) ────────────────
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf "╔══════════════════════════════════════════════════════════════╗\n"
  printf "║  mammalia_nt_db.sh — SLURM submission                       ║\n"
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

  # Prompt 2: download/refresh nt_euk
  printf "\nnt_euk source database options:\n"
  printf "  [1] Use existing nt_euk (skip download) — recommended if already present\n"
  printf "  [2] Download / refresh nt_euk from NCBI via rsync (~400 GB, several hours)\n"
  read -rp "Choice [1/2, default=1]: " DL_CHOICE
  case "${DL_CHOICE:-1}" in
    2) DOWNLOAD_DB=true  ;;
    *) DOWNLOAD_DB=false ;;
  esac
  printf "[INFO] DOWNLOAD_DB=%s\n" "$DOWNLOAD_DB"

  mkdir -p logs

  JOB_ID=$(DOWNLOAD_DB="$DOWNLOAD_DB" \
    sbatch --account="$ACCOUNT" "$(realpath "$0")" | awk '{print $NF}')

  printf "[OK]   Submitted job %s with account %s\n" "$JOB_ID" "$ACCOUNT"
  printf "       Monitor : squeue -j %s\n" "$JOB_ID"
  printf "       Log     : logs/mammalia_nt_%s.out\n" "$JOB_ID"
  exit 0
fi

# ─── Everything below runs inside the SLURM job ──────────────────────────────

set -eo pipefail
IFS=$'\n\t'

# TACC's tacc-apptainer bash-completion script references BASH_COMPLETION_DEBUG
# without a default; under `set -u` this aborts the job during module load.
# Enable nounset after modules are loaded.
module load tacc-apptainer
module load biocontainers
module load ncbi-blast-plus

set -u

export LC_ALL=C

BLASTDB_DIR="$WORK/Metabarcoding/blastdb"
DB_PREFIX="${BLASTDB_DIR}/nt_euk"
OUT_PREFIX="${BLASTDB_DIR}/mammalia_nt"
STALE_DAYS="${STALE_DAYS:-90}"              # age threshold for staleness warning
DOWNLOAD_DB="${DOWNLOAD_DB:-false}"         # must be explicitly set to "true"

mkdir -p "$BLASTDB_DIR" logs

printf "╔══════════════════════════════════════════════════════════════╗\n"
printf "║  mammalia_nt_db.sh — Mammalia BLAST DB builder              ║\n"
printf "╚══════════════════════════════════════════════════════════════╝\n"
printf "[INFO] BLASTDB dir  : %s\n" "$BLASTDB_DIR"
printf "[INFO] DOWNLOAD_DB  : %s\n" "$DOWNLOAD_DB"

export BLASTDB="$BLASTDB_DIR"

# ─── Helper: find the newest nt_euk .nhr index file ──────────────────────────
_newest_nhr() {
  local newest=0 ts
  for f in "${DB_PREFIX}".*.nhr "${DB_PREFIX}".nhr; do
    [[ -f "$f" ]] || continue
    ts=$(stat -c %Y "$f" 2>/dev/null || echo 0)
    (( ts > newest )) && newest=$ts
  done
  echo "$newest"
}

# ─── Helper: check whether nt_euk index files exist ──────────────────────────
_nt_euk_present() {
  ls "${DB_PREFIX}".*.nhr &>/dev/null || ls "${DB_PREFIX}".nhr &>/dev/null
}

# ─── Download function (only called when DOWNLOAD_DB=true) ───────────────────
_download_nt_euk() {
  printf "[INFO] Downloading nt_euk volumes from NCBI via rsync...\n"
  cd "$BLASTDB_DIR"

  rsync --list-only --no-motd rsync://ftp.ncbi.nlm.nih.gov/blast/db/ \
    | awk '/nt_euk.*\.tar\.gz$/ {print $NF}' > nt_euk.links.list

  local n_files
  n_files=$(wc -l < nt_euk.links.list)
  printf "[INFO] Found %s volume files to download\n" "$n_files"

  rsync --archive --partial --progress --compress --no-motd \
    --files-from=nt_euk.links.list \
    rsync://ftp.ncbi.nlm.nih.gov/blast/db/ .

  printf "[INFO] Decompressing volumes...\n"
  for f in nt_euk.*.tar.gz; do
    [[ -f "$f" ]] || continue
    tar -zxf "$f" && rm -f "$f"
  done

  printf "[OK]   nt_euk download and decompression complete\n"
}

# ═══════════════════════════════════════════════════════════════════════════════
# Step 1: Resolve nt_euk source database
# ═══════════════════════════════════════════════════════════════════════════════

if [[ "$DOWNLOAD_DB" == "true" ]]; then
  # ── User explicitly opted in to download ──
  printf "[INFO] DOWNLOAD_DB=true — downloading/refreshing nt_euk\n"
  _download_nt_euk

elif _nt_euk_present; then
  # ── nt_euk exists — use it, but warn if stale ──
  NEWEST=$(_newest_nhr)
  NOW=$(date +%s)
  AGE_DAYS=$(( (NOW - NEWEST) / 86400 ))

  if (( AGE_DAYS > STALE_DAYS )); then
    printf "[WARN] nt_euk index files are %s days old (> %s days)\n" \
      "$AGE_DAYS" "$STALE_DAYS" >&2
    printf "       To refresh, resubmit with:\n" >&2
    printf "       DOWNLOAD_DB=true sbatch scripts/mammalia_nt_db.sh\n" >&2
  fi
  printf "[OK]   Using existing nt_euk (%s days old)\n" "$AGE_DAYS"

else
  # ── nt_euk missing and user did NOT opt in — stop ──
  printf "[ERROR] nt_euk index files not found in %s\n" "$BLASTDB_DIR" >&2
  printf "        The source database must be downloaded before building\n" >&2
  printf "        the Mammalia DB. Resubmit with:\n" >&2
  printf "\n" >&2
  printf "          DOWNLOAD_DB=true sbatch scripts/mammalia_nt_db.sh\n" >&2
  printf "\n" >&2
  exit 1
fi

# ─── Verify nt_euk is usable ─────────────────────────────────────────────────
if ! _nt_euk_present; then
  printf "[ERROR] nt_euk index files not found after download step\n" >&2
  exit 1
fi
printf "[OK]   nt_euk index files verified\n"

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Extract Mammalia sequences, excluding Hominidae
# ═══════════════════════════════════════════════════════════════════════════════

printf "[INFO] Extracting Mammalia (taxid 40674) excluding Hominidae (taxid 9604)...\n"
blastdbcmd \
  -db "${DB_PREFIX}" \
  -taxids 40674 \
  -negative_taxids 9604 \
  -outfmt "%f" \
  -out "${BLASTDB_DIR}/nt_euk_mammalia.fasta"

NSEQ=$(grep -c '^>' "${BLASTDB_DIR}/nt_euk_mammalia.fasta" 2>/dev/null || echo 0)
printf "[OK]   Extracted %s sequences\n" "$NSEQ"

if [[ "$NSEQ" -eq 0 ]]; then
  printf "[ERROR] Extraction produced 0 sequences — check nt_euk integrity\n" >&2
  exit 1
fi

# ═══════════════════════════════════════════════════════════════════════════════
# Step 3: Build Mammalia-only BLAST database
# ═══════════════════════════════════════════════════════════════════════════════

printf "[INFO] Building BLAST DB: mammalia_nt...\n"
makeblastdb \
  -in "${BLASTDB_DIR}/nt_euk_mammalia.fasta" \
  -dbtype nucl \
  -parse_seqids \
  -out "${OUT_PREFIX}" \
  -title "Mammalia_nt_excluding_Hominidae"

printf "[OK]   Mammalia NT database ready: %s.*\n" "$OUT_PREFIX"

# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Cleanup intermediate FASTA
# ═══════════════════════════════════════════════════════════════════════════════

rm -f "${BLASTDB_DIR}/nt_euk_mammalia.fasta"
printf "[OK]   Intermediate FASTA removed\n"
printf "[OK]   Done. Run: bash scripts/setup.sh eDNA/2024 blastdb/mammalia_nt 8 \"\" otu\n"
printf "       Then submit: sbatch scripts/batch_job.sh\n"
