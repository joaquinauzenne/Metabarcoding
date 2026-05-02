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
#
# Note: we deliberately do NOT load `biocontainers` / `ncbi-blast-plus`.
# The system's `/usr/bin/blastdbcmd` resolved via those modules is XALT-
# wrapped and links against a newer glibc than the LS6 compute node
# provides (`GLIBC_2.25' not found`). Instead we run BLAST+ inside the
# edna-pipeline.sif container, which is the pattern the rest of the
# pipeline already uses.
module load tacc-apptainer

set -u

SIF="$WORK/Metabarcoding/Apptainer/edna-pipeline.sif"
if [[ ! -f "$SIF" ]]; then
  printf "[ERROR] edna-pipeline.sif not found at %s\n" "$SIF" >&2
  exit 1
fi

export LC_ALL=C

BLASTDB_DIR="$WORK/Metabarcoding/blastdb"
DB_PREFIX="${BLASTDB_DIR}/nt_euk"
OUT_PREFIX="${BLASTDB_DIR}/mammalia_nt"
STALE_DAYS="${STALE_DAYS:-90}"              # age threshold for staleness warning

# ─── Prompt for DOWNLOAD_DB if unset (only in interactive mode) ────────────────
if [[ -z "${DOWNLOAD_DB:-}" ]]; then
  if [[ -t 0 ]]; then
    # Interactive: stdin is a TTY, prompt the user
    printf "\nDownload / refresh nt_euk from NCBI? (~400 GB, several hours)\n"
    printf "  [1] Use existing nt_euk (skip download)\n"
    printf "  [2] Download / refresh nt_euk from NCBI via rsync\n"
    read -rp "Choice [1/2, default=1]: " DL_CHOICE
    case "${DL_CHOICE:-1}" in
      2) DOWNLOAD_DB=true  ;;
      *) DOWNLOAD_DB=false ;;
    esac
  else
    # Non-interactive (batch job): stdin is not a TTY, require explicit env var
    printf "[ERROR] DOWNLOAD_DB not set and stdin is not a TTY (batch mode).\n" >&2
    printf "        Set DOWNLOAD_DB before submitting:\n" >&2
    printf "          DOWNLOAD_DB=true sbatch scripts/mammalia_nt_db.sh\n" >&2
    printf "        Or run interactively in idev:\n" >&2
    printf "          idev -t 00:15:00 -p development\n" >&2
    printf "          bash scripts/mammalia_nt_db.sh\n" >&2
    exit 1
  fi
else
  DOWNLOAD_DB="${DOWNLOAD_DB}"              # already set, use as-is
fi

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

# ─── Helper: check whether nt_euk volume index files exist ───────────────────
_nt_euk_present() {
  ls "${DB_PREFIX}".*.nhr &>/dev/null || ls "${DB_PREFIX}".nhr &>/dev/null
}

# ─── Helper: list missing auxiliary files required by -taxids filtering ──────
# blastdbcmd -taxids needs:
#   nt_euk.nal              — multi-volume alias (without it, blastdbcmd can't
#                             resolve "-db nt_euk" even when volumes exist)
#   taxdb.btd / taxdb.bti   — scientific-name / common-name lookup tables
#   taxonomy4blast.sqlite3  — lineage DB for taxid subtree expansion (required
#                             to resolve "all descendants of taxid 40674")
_missing_aux_files() {
  local missing=()
  [[ -f "${BLASTDB_DIR}/nt_euk.nal" ]]             || missing+=("nt_euk.nal")
  [[ -f "${BLASTDB_DIR}/taxdb.btd" ]]              || missing+=("taxdb.btd")
  [[ -f "${BLASTDB_DIR}/taxdb.bti" ]]              || missing+=("taxdb.bti")
  [[ -f "${BLASTDB_DIR}/taxonomy4blast.sqlite3" ]] || missing+=("taxonomy4blast.sqlite3")
  if (( ${#missing[@]} > 0 )); then
    printf "%s\n" "${missing[@]}"
  fi
}

# ─── Download function (only called when DOWNLOAD_DB=true) ───────────────────
# Fetches four classes of files via rsync, then decompresses:
#   1. nt_euk.NN.tar.gz          — per-volume index archives
#   2. nt_euk.nal                — multi-volume alias (standalone on NCBI FTP)
#   3. taxdb.tar.gz              — expands to taxdb.btd + taxdb.bti
#   4. taxonomy4blast.sqlite3.gz — expands to taxonomy4blast.sqlite3
_download_nt_euk() {
  printf "[INFO] Downloading nt_euk + taxonomy files from NCBI via rsync...\n"
  cd "$BLASTDB_DIR"

  rsync --list-only --no-motd rsync://ftp.ncbi.nlm.nih.gov/blast/db/ \
    | awk '
        $NF ~ /^nt_euk\.[0-9]+\.tar\.gz$/ { print $NF }
        $NF == "nt_euk.nal"                { print $NF }
        $NF == "taxdb.tar.gz"              { print $NF }
        $NF == "taxonomy4blast.sqlite3.gz" { print $NF }
      ' > nt_euk.links.list

  local n_files
  n_files=$(wc -l < nt_euk.links.list)
  printf "[INFO] Found %s file(s) to download (volumes + alias + taxonomy)\n" "$n_files"

  rsync --archive --partial --progress --compress --no-motd \
    --files-from=nt_euk.links.list \
    rsync://ftp.ncbi.nlm.nih.gov/blast/db/ .

  printf "[INFO] Decompressing nt_euk volume archives...\n"
  for f in nt_euk.*.tar.gz; do
    [[ -f "$f" ]] || continue
    tar -zxf "$f" && rm -f "$f"
  done

  if [[ -f taxdb.tar.gz ]]; then
    printf "[INFO] Decompressing taxdb.tar.gz...\n"
    tar -zxf taxdb.tar.gz && rm -f taxdb.tar.gz
  fi

  if [[ -f taxonomy4blast.sqlite3.gz ]]; then
    printf "[INFO] Decompressing taxonomy4blast.sqlite3.gz...\n"
    gunzip -f taxonomy4blast.sqlite3.gz
  fi

  printf "[OK]   nt_euk + taxonomy download and decompression complete\n"
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
  printf "[ERROR] nt_euk volume files (.nhr) not found after download step\n" >&2
  exit 1
fi

# blastdbcmd -taxids 40674 will fail without the alias + taxonomy support
# files (seen empirically: "Unknown database" or "-taxids requires additional
# data files"). Check them up front so we fail fast with actionable guidance.
MISSING_AUX=$(_missing_aux_files)
if [[ -n "$MISSING_AUX" ]]; then
  printf "[ERROR] Missing required BLAST auxiliary files in %s:\n" "$BLASTDB_DIR" >&2
  while IFS= read -r f; do
    printf "          - %s\n" "$f" >&2
  done <<< "$MISSING_AUX"
  printf "\n" >&2
  printf "        blastdbcmd -taxids needs:\n" >&2
  printf "          nt_euk.nal             (multi-volume alias)\n" >&2
  printf "          taxdb.btd / taxdb.bti  (name lookup)\n" >&2
  printf "          taxonomy4blast.sqlite3 (lineage expansion)\n" >&2
  printf "\n" >&2
  printf "        To fetch missing files, resubmit with:\n\n" >&2
  printf "          DOWNLOAD_DB=true sbatch scripts/mammalia_nt_db.sh\n\n" >&2
  exit 1
fi
printf "[OK]   nt_euk volumes, alias, and taxonomy files verified\n"

# ═══════════════════════════════════════════════════════════════════════════════
# Step 2: Extract Mammalia sequences (Hominidae excluded at search time)
# ═══════════════════════════════════════════════════════════════════════════════
# blastdbcmd supports only -taxids/-taxidlist (no negative variant). Hominidae
# (9604) is excluded downstream by blastn -negative_taxids in pipeline.sh's
# FILTER_OPT (Phase 4), so this DB intentionally contains all of Mammalia.

printf "[INFO] Extracting Mammalia (taxid 40674); Hominidae (9604) excluded at search time via blastn -negative_taxids\n"
apptainer exec \
  --bind "$WORK":"$WORK" \
  --env "BLASTDB=$BLASTDB_DIR" \
  "$SIF" \
  blastdbcmd \
    -db "${DB_PREFIX}" \
    -taxids 40674 \
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
apptainer exec \
  --bind "$WORK":"$WORK" \
  "$SIF" \
  makeblastdb \
    -in "${BLASTDB_DIR}/nt_euk_mammalia.fasta" \
    -dbtype nucl \
    -parse_seqids \
    -out "${OUT_PREFIX}" \
    -title "Mammalia_nt (Hominidae excluded at query time)"

printf "[OK]   Mammalia NT database ready: %s.*\n" "$OUT_PREFIX"

# ═══════════════════════════════════════════════════════════════════════════════
# Step 4: Cleanup intermediate FASTA
# ═══════════════════════════════════════════════════════════════════════════════

rm -f "${BLASTDB_DIR}/nt_euk_mammalia.fasta"
printf "[OK]   Intermediate FASTA removed\n"
printf "[OK]   Done. Run: bash scripts/setup.sh eDNA/2024 blastdb/mammalia_nt 8 \"\" otu\n"
printf "       Then submit: sbatch scripts/batch_job.sh\n"
