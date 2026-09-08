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

# ── Load pipeline configuration (same source of truth as pipeline.sh) ────────
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_CONFIG="${PIPELINE_CONFIG:-$SCRIPTS_DIR/../config/pipeline.yaml}"
if [[ ! -f "$PIPELINE_CONFIG" ]]; then
  printf "[ERROR] Pipeline config not found: %s\n" "$PIPELINE_CONFIG" >&2
  printf "        Set PIPELINE_CONFIG or create config/pipeline.yaml\n" >&2
  exit 1
fi
# shellcheck source=lib/parse_yaml.sh
source "$SCRIPTS_DIR/lib/parse_yaml.sh"
eval "$(parse_yaml "$PIPELINE_CONFIG" CFG_)"

printf "[INFO] eDNA pipeline preflight\n"
printf "  RAW_DIR:  %s\n" "$RAW_DIR"
printf "  DB_DIR:   %s\n" "$DB_DIR"
printf "  THREADS:  %s\n" "$THREADS"
printf "  MODE:     %s\n" "$MODE"
printf "  CONFIG:   %s\n" "$PIPELINE_CONFIG"

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
# DB_DIR is a DIRECTORY holding the database; the DB stem comes from the config
# (blast.db_name). This mirrors pipeline.sh, which queries "$DB_DIR/$db_name".
# Checking DB_DIR as if it were a path prefix (the old behaviour) reported a
# spurious WARN even when a valid database was present.
# ─────────────────────────────────────────────────────────────────────────────
DB_NAME="${CFG_blast_db_name:-}"
if [[ -z "$DB_NAME" ]]; then
  printf "[ERROR] blast.db_name missing from %s\n" "$PIPELINE_CONFIG" >&2; exit 1
fi
DB_PREFIX="$DB_DIR/$DB_NAME"
printf "[INFO] Checking BLAST database: %s (dir=%s, db_name=%s)\n" "$DB_PREFIX" "$DB_DIR" "$DB_NAME"

if [[ ! -d "$DB_DIR" ]]; then
  printf "[ERROR] DB directory '%s' does not exist\n" "$DB_DIR" >&2; exit 1
fi
if compgen -G "${DB_PREFIX}.*n[hsi]r" >/dev/null || compgen -G "${DB_PREFIX}*.nal" >/dev/null; then
  printf "[OK]   BLAST DB found: %s.*\n" "$DB_PREFIX"
else
  printf "[WARN] No BLAST DB index files found for '%s.*'\n" "$DB_PREFIX" >&2
  printf "       Databases present in %s:\n" "$DB_DIR" >&2
  { compgen -G "$DB_DIR/*.nal"; compgen -G "$DB_DIR/*.nin"; } 2>/dev/null \
    | sed -e 's/\.nal$//' -e 's/\.nin$//' -e 's/\.[0-9][0-9]*$//' \
    | xargs -r -I{} basename {} | sort -u | sed 's/^/         /' >&2 || true
  printf "       Build one with scripts/build_midori_db.sh (MIDORI2).\n" >&2
fi

# Taxonomy lookup files: required for -taxids filtering and for populating the
# staxids/sscinames columns Phase 4 and lca_assignment.R depend on.
for tf in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do
  if [[ -f "$DB_DIR/$tf" ]]; then
    printf "[OK]   taxonomy file: %s\n" "$tf"
  else
    printf "[WARN] Missing %s in %s -- '-taxids' filtering and sscinames will fail\n" "$tf" "$DB_DIR" >&2
  fi
done

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
# 7. batch_job.sh allocation handling
# batch_job.sh no longer hardcodes an allocation — it prompts for the TACC
# account interactively at submit time (read -rp) and self-submits via sbatch.
# Confirm that prompt is present rather than grepping for a static '#SBATCH -A'
# line: that grep matches nothing now and, piped under `set -euo pipefail`, would
# abort setup.sh before the manifest step below.
# ─────────────────────────────────────────────────────────────────────────────
if grep -q 'Enter your TACC allocation account' scripts/batch_job.sh 2>/dev/null; then
  printf "[OK]   batch_job.sh prompts for the TACC allocation at submit time\n"
elif grep -q '#SBATCH -A' scripts/batch_job.sh 2>/dev/null; then
  ALLOC="$(grep '#SBATCH -A' scripts/batch_job.sh | awk '{print $3}')"
  printf "[OK]   batch_job.sh allocation: %s\n" "$ALLOC"
else
  printf "[WARN] batch_job.sh has no interactive allocation prompt and no '#SBATCH -A' — verify before sbatch\n" >&2
fi

# ─────────────────────────────────────────────────────────────────────────────
# 8. Provenance manifest (reproducibility)
# Emits logs/run_manifest_<ts>.json + logs/env_snapshot_<ts>.txt capturing the
# exact state of this run: git commit, config hash, container env, and input
# FASTQ inventory. Best-effort — individual failures warn but do not abort.
# Set MANIFEST_CHECKSUMS=false to skip per-file sha256 (faster on large inputs).
# ─────────────────────────────────────────────────────────────────────────────
printf "[INFO] Writing provenance manifest...\n"
TS="$(date -u +%Y%m%dT%H%M%SZ)"
MANIFEST="logs/run_manifest_${TS}.json"
ENV_SNAPSHOT="logs/env_snapshot_${TS}.txt"
CONFIG_FILE="$PIPELINE_CONFIG"

# Container environment snapshot — the true pinned versions from the built .sif.
if apptainer exec "$SIF" conda list -n edna-pipeline --export > "$ENV_SNAPSHOT" 2>/dev/null; then
  printf "[OK]   Env snapshot: %s\n" "$ENV_SNAPSHOT"
else
  printf "[WARN] Could not capture conda env snapshot from container\n" >&2
  ENV_SNAPSHOT=""
fi

GIT_COMMIT="$(git rev-parse HEAD 2>/dev/null || echo unknown)"
if git diff --quiet 2>/dev/null && git diff --cached --quiet 2>/dev/null; then
  GIT_DIRTY=false
else
  GIT_DIRTY=true
fi

if [[ -f "$CONFIG_FILE" ]]; then
  CONFIG_SHA="$(sha256sum "$CONFIG_FILE" 2>/dev/null | awk '{print $1}')"
else
  CONFIG_SHA=""
fi

HOSTN="$(hostname)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then EXEC_CTX="compute"
elif [[ "$HOSTN" == login* ]]; then EXEC_CTX="login"
else EXEC_CTX="local"; fi

DO_SUM="${MANIFEST_CHECKSUMS:-true}"
{
  printf '{\n'
  printf '  "generated_at": "%s",\n' "$TS"
  printf '  "host": "%s",\n' "$HOSTN"
  printf '  "execution_context": "%s",\n' "$EXEC_CTX"
  printf '  "slurm_job_id": "%s",\n' "${SLURM_JOB_ID:-}"
  printf '  "git_commit": "%s",\n' "$GIT_COMMIT"
  printf '  "git_dirty": %s,\n' "$GIT_DIRTY"
  printf '  "config_path": "%s",\n' "$CONFIG_FILE"
  printf '  "config_sha256": "%s",\n' "$CONFIG_SHA"
  printf '  "container_sif": "%s",\n' "$SIF"
  printf '  "env_snapshot": "%s",\n' "$ENV_SNAPSHOT"
  printf '  "mode": "%s",\n' "$MODE"
  printf '  "threads": "%s",\n' "$THREADS"
  printf '  "db_dir": "%s",\n' "$DB_DIR"
  printf '  "raw_dir": "%s",\n' "$RAW_DIR"
  printf '  "raw_fastq": [\n'
  mapfile -t FQ < <(find "$RAW_DIR" -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) | sort)
  for i in "${!FQ[@]}"; do
    f="${FQ[$i]}"
    sz="$(stat -c %s "$f" 2>/dev/null || echo 0)"
    if [[ "$DO_SUM" == "true" ]]; then
      sh="$(sha256sum "$f" 2>/dev/null | awk '{print $1}')"
    else
      sh=""
    fi
    sep=","
    [[ "$i" -eq $(( ${#FQ[@]} - 1 )) ]] && sep=""
    printf '    {"file": "%s", "bytes": %s, "sha256": "%s"}%s\n' "$f" "$sz" "$sh" "$sep"
  done
  printf '  ]\n'
  printf '}\n'
} > "$MANIFEST"
printf "[OK]   Manifest: %s (%d input FASTQ, checksums=%s)\n" "$MANIFEST" "${#FQ[@]}" "$DO_SUM"

printf "[OK]   Setup complete.\n"
