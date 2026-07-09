#!/usr/bin/env bash
set -eo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# scripts/build_midori_db.sh  — M2 (Phase B1): build the MIDORI2 reference DB
#
# Purpose: download MIDORI2 UNIQ (BLAST-flavour) reference FASTAs for 12S (srRNA)
#          and 16S (lrRNA), build a single BLAST database with makeblastdb
#          -parse_seqids, and record exact provenance. Reuses the download +
#          makeblastdb path already PROVEN by scripts/diag_midori_blast.sh
#          (makeblastdb -parse_seqids + the co-located taxdb files resolve
#          staxids/sscinames at query time).
#
# Mammalia scoping: this DB contains the MIDORI2 12S+16S reference as-is. The
# pipeline restricts hits to Mammalia at QUERY time via `-taxids 40674
# -negative_taxids 9604` (see pipeline.sh Phase 4; proven with the taxdb files).
# We do NOT filter by parsing FASTA headers (their exact lineage format is not
# assumed here). An optional physical subset via `blastdbcmd -taxidlist` is left
# as a documented enhancement — the query-time filter already yields mammal-only
# results without it.
#
# Capability / limitation (paired): DOES produce a versioned, citable mito
# reference whose accessions resolve to taxa via taxdb. Does NOT expand Neotropical
# coverage — run db_coverage_audit.R to quantify how much of your checklist the
# release actually covers.
#
# Run in an idev session (download = network):
#   idev -t 00:30:00 -p development
#   module load tacc-apptainer
#   MIDORI_12S_URL="<srRNA UNIQ BLAST fasta url>" \
#   MIDORI_16S_URL="<lrRNA UNIQ BLAST fasta url>" \
#   RELEASE="GB270" \
#   bash scripts/build_midori_db.sh
#
# Get current-release URLs from https://www.reference-midori.info/download.php
# (UNIQ flavour, BLAST format; srRNA = 12S, lrRNA = 16S). URLs may be .fasta.gz
# or .fasta.zip — both are handled.
#
# Env:
#   MIDORI_12S_URL   (required) 12S/srRNA UNIQ BLAST FASTA url
#   MIDORI_16S_URL   (required) 16S/lrRNA UNIQ BLAST FASTA url
#   OUT_DIR          DB output dir [default: $WORK/Metabarcoding/blastdb]
#                    (default co-locates with the retained taxdb.* files)
#   DB_NAME          makeblastdb -out stem [default: mammalia_midori]
#                    (MUST match blast.db_name in config/pipeline.yaml)
#   RELEASE          release label for provenance [default: unknown]
#   SIF              container [default: Apptainer/edna-pipeline.sif]
# ─────────────────────────────────────────────────────────────────────────────

MIDORI_12S_URL="${MIDORI_12S_URL:-}"
MIDORI_16S_URL="${MIDORI_16S_URL:-}"
OUT_DIR="${OUT_DIR:-${WORK:-$PWD}/Metabarcoding/blastdb}"
DB_NAME="${DB_NAME:-mammalia_midori}"
RELEASE="${RELEASE:-unknown}"
SIF="${SIF:-Apptainer/edna-pipeline.sif}"
BIND_PATH="${WORK:-$PWD}"

WORK_TMP="$OUT_DIR/midori_build"
COMBINED="$WORK_TMP/midori_12S16S.fasta"
DB_PREFIX="$OUT_DIR/$DB_NAME"
BUILD_JSON="${DB_PREFIX}.build.json"

# ─── Guards ──────────────────────────────────────────────────────────────────
HOSTNAME_VAL="$(hostname)"
if [[ "$HOSTNAME_VAL" == login* ]] && [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf "[ERROR] Refusing to run on login node '%s'. Use idev.\n" "$HOSTNAME_VAL" >&2; exit 1
fi
[[ -n "$MIDORI_12S_URL" ]] || { printf "[ERROR] MIDORI_12S_URL not set (see header for where to find it)\n" >&2; exit 1; }
[[ -n "$MIDORI_16S_URL" ]] || { printf "[ERROR] MIDORI_16S_URL not set\n" >&2; exit 1; }
command -v apptainer &>/dev/null || { printf "[ERROR] apptainer not found. module load tacc-apptainer\n" >&2; exit 1; }
[[ -f "$SIF" ]] || { printf "[ERROR] Container not found: %s\n" "$SIF" >&2; exit 1; }

mkdir -p "$WORK_TMP"

# ─── Download + decompress helper (reuses diag_midori_blast.sh proven logic) ──
# Handles .fasta.gz and .fasta.zip; writes a plain .fasta to $2.
download_fasta() {
  local url="$1" out_fasta="$2" dl fmt
  case "$url" in
    *.fasta.gz)  dl="${out_fasta}.gz";  fmt="gz"  ;;
    *.fasta.zip) dl="${out_fasta}.zip"; fmt="zip" ;;
    *.zip)       dl="${out_fasta}.zip"; fmt="zip" ;;
    *)           dl="$out_fasta";       fmt="none" ;;
  esac
  if [[ -s "$out_fasta" ]]; then
    printf "[SKIP] Already present: %s\n" "$out_fasta"; return 0
  fi
  printf "[INFO] Downloading: %s\n" "$url"
  curl -fSL --retry 3 --retry-delay 5 -o "$dl" "$url"
  case "$fmt" in
    gz)  gunzip -c "$dl" > "$out_fasta" ;;
    zip) # Extract into an isolated dir so a previously-built marker FASTA already
         # sitting in $WORK_TMP (e.g. midori_12S.fasta) can never be picked up here.
         local exdir extracted
         exdir="$(mktemp -d "${WORK_TMP}/unzip.XXXXXX")"
         unzip -o -j "$dl" "*.fasta" -d "$exdir" >/dev/null
         extracted="$(find "$exdir" -maxdepth 1 -name '*.fasta' | head -1 || true)"
         [[ -n "$extracted" ]] && mv "$extracted" "$out_fasta"
         rm -rf "$exdir" ;;
    none) : ;;
  esac
  [[ -s "$out_fasta" ]] || { printf "[ERROR] Download/decompress produced no FASTA: %s\n" "$out_fasta" >&2; return 1; }
}

FASTA_12S="$WORK_TMP/midori_12S.fasta"
FASTA_16S="$WORK_TMP/midori_16S.fasta"

printf "[INFO] === Step 1: download MIDORI2 12S + 16S references ===\n"
download_fasta "$MIDORI_12S_URL" "$FASTA_12S"
download_fasta "$MIDORI_16S_URL" "$FASTA_16S"

N12=$(grep -c '^>' "$FASTA_12S" || echo 0)
N16=$(grep -c '^>' "$FASTA_16S" || echo 0)
printf "[INFO] 12S sequences: %s ; 16S sequences: %s\n" "$N12" "$N16"
(( N12 > 0 && N16 > 0 )) || { printf "[ERROR] One of the marker FASTAs is empty\n" >&2; exit 1; }

printf "[INFO] Sample 12S header: %s\n" "$(grep -m1 '^>' "$FASTA_12S" || true)"
printf "[INFO] Sample 16S header: %s\n" "$(grep -m1 '^>' "$FASTA_16S" || true)"

# ─── Step 2: concatenate ─────────────────────────────────────────────────────
printf "\n[INFO] === Step 2: concatenate 12S + 16S -> %s ===\n" "$COMBINED"
cat "$FASTA_12S" "$FASTA_16S" > "$COMBINED"
NTOT=$(grep -c '^>' "$COMBINED" || echo 0)
printf "[INFO] Combined sequences: %s\n" "$NTOT"

# ─── Step 3: makeblastdb -parse_seqids (proven path) ─────────────────────────
printf "\n[INFO] === Step 3: makeblastdb -> %s ===\n" "$DB_PREFIX"
apptainer exec --bind "$BIND_PATH":"$BIND_PATH" "$SIF" \
  makeblastdb \
    -in "$COMBINED" \
    -dbtype nucl \
    -parse_seqids \
    -out "$DB_PREFIX" \
    -title "MIDORI2 12S+16S ${RELEASE} (mammal ID via query-time -taxids 40674)"
printf "[OK]   makeblastdb complete: %s.*\n" "$DB_PREFIX"

# ─── Step 4: self-query verification (staxids/sscinames populate?) ───────────
# Same success criterion diag_midori_blast.sh proved: pick the first sequence,
# BLAST it against the new DB, confirm taxonomy columns are non-empty.
printf "\n[INFO] === Step 4: verify taxonomy resolution (self-query) ===\n"
TEST_QUERY="$WORK_TMP/test_query.fasta"
awk '/^>/{if(n++)exit; print; next} n' "$COMBINED" > "$TEST_QUERY"   # file-read, no SIGPIPE
BLAST_OUT="$WORK_TMP/verify.tsv"
apptainer exec --bind "$BIND_PATH":"$BIND_PATH" --env "BLASTDB=$OUT_DIR" "$SIF" \
  blastn -query "$TEST_QUERY" -db "$DB_PREFIX" \
    -outfmt "6 qseqid sseqid pident staxids sscinames" \
    -max_target_seqs 3 -perc_identity 90 -out "$BLAST_OUT" 2>/dev/null || true

VERIFY_STATUS="FAIL"
if [[ -s "$BLAST_OUT" ]]; then
  EMPTY_TAX=$(awk -F'\t' '$4=="" || $4=="N/A"' "$BLAST_OUT" | wc -l)
  if (( EMPTY_TAX == 0 )); then
    VERIFY_STATUS="PASS"
    printf "[OK]   Taxonomy columns populated. First hit: %s\n" "$(head -1 "$BLAST_OUT")"
  else
    printf "[WARN] %d hit(s) have empty staxids — check that taxdb.* live in %s\n" "$EMPTY_TAX" "$OUT_DIR" >&2
  fi
else
  printf "[WARN] Self-query returned no hits — verification inconclusive\n" >&2
fi

# ─── Step 5: provenance record ───────────────────────────────────────────────
sha() { sha256sum "$1" 2>/dev/null | awk '{print $1}'; }
BLAST_VER="$(apptainer exec "$SIF" blastn -version 2>/dev/null | head -1 || echo unknown)"
{
  printf '{\n'
  printf '  "db_name": "%s",\n' "$DB_NAME"
  printf '  "release": "%s",\n' "$RELEASE"
  printf '  "built_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '  "blast_version": "%s",\n' "$BLAST_VER"
  printf '  "midori_12S_url": "%s",\n' "$MIDORI_12S_URL"
  printf '  "midori_16S_url": "%s",\n' "$MIDORI_16S_URL"
  printf '  "n_seqs_12S": %s,\n' "$N12"
  printf '  "n_seqs_16S": %s,\n' "$N16"
  printf '  "n_seqs_total": %s,\n' "$NTOT"
  printf '  "sha256_12S": "%s",\n' "$(sha "$FASTA_12S")"
  printf '  "sha256_16S": "%s",\n' "$(sha "$FASTA_16S")"
  printf '  "mammalia_scoping": "query-time -taxids 40674 -negative_taxids 9604 (pipeline.sh)",\n'
  printf '  "taxonomy_self_query": "%s"\n' "$VERIFY_STATUS"
  printf '}\n'
} > "$BUILD_JSON"

printf "\n[OK]   MIDORI2 DB built: %s.*\n" "$DB_PREFIX"
printf "       Provenance:       %s\n" "$BUILD_JSON"
printf "       Verification:     %s\n" "$VERIFY_STATUS"
printf "       Point pipeline.sh at DB_DIR=%s (blast.db_name=%s in config).\n" "$OUT_DIR" "$DB_NAME"
printf "       Next: scripts/db_coverage_audit.R to quantify checklist coverage.\n"
[[ "$VERIFY_STATUS" == "PASS" ]] || printf "[WARN] Taxonomy verification did not PASS — resolve before trusting Phase 4.\n" >&2
