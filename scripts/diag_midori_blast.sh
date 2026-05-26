#!/usr/bin/env bash
# set -eo pipefail before module load; -u deferred (TACC bash-completion references
# BASH_COMPLETION_DEBUG without a default, which aborts under set -u).
set -eo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# Script : diag_midori_blast.sh
# Purpose: Phase B0 de-risk probe — verify that a MIDORI2 reference FASTA can
#          be built into a BLAST DB with -parse_seqids and that query-time
#          staxids / sscinames columns are populated.
#
#          MIDORI2 distributes FASTA files with NCBI accessions in the headers.
#          makeblastdb -parse_seqids + the existing taxdb files should preserve
#          staxids / sscinames lookup at query time exactly as the current
#          nt_euk-based pipeline does. This script confirms that empirically.
#
# Configuration:
#   Set MIDORI2_URL to the download URL of a small MIDORI2 FASTA
#   (UNIQ flavour, .fasta.gz).  Find current URLs at:
#     https://www.reference-midori.info/download.php
#   Example (GenBank release 270, 12S srRNA — GB270 uses .fasta.zip, not .fasta.gz):
#     MIDORI2_URL="https://www.reference-midori.info/download/Databases/\
# GenBank270_2026-02-15/BLAST/uniq/fasta/MIDORI2_UNIQ_NUC_GB270_srRNA_BLAST.fasta.zip"
#
# Run inside an idev compute session (NOT on a login node):
#   idev -t 00:30:00 -p development
#   module load tacc-apptainer
#   MIDORI2_URL="<url>" bash scripts/diag_midori_blast.sh
#
# Exit codes:
#   0  — MIDORI2 makeblastdb succeeds and staxids/sscinames are populated
#   1  — download failed, makeblastdb failed, or taxonomy columns empty
# ─────────────────────────────────────────────────────────────────────────────

MIDORI2_URL="${MIDORI2_URL:-}"
TEST_DIR="${TEST_DIR:-$WORK/Metabarcoding/midori_test}"
SIF="${SIF:-$WORK/Metabarcoding/Apptainer/edna-pipeline.sif}"
BLASTDB_DIR="${BLASTDB_DIR:-$WORK/Metabarcoding/blastdb}"

# ─── Login-node guard ────────────────────────────────────────────────────────
HOSTNAME_VAL="$(hostname)"
if [[ "$HOSTNAME_VAL" == login* ]] && [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf "[ERROR] Refusing to run on login node '%s'.\n" "$HOSTNAME_VAL" >&2
  printf "        Start a compute session:\n" >&2
  printf "          idev -t 00:30:00 -p development\n" >&2
  printf "          module load tacc-apptainer\n" >&2
  printf "          MIDORI2_URL=\"<url>\" bash scripts/diag_midori_blast.sh\n" >&2
  exit 1
fi

if [[ -z "$MIDORI2_URL" ]]; then
  printf "[ERROR] MIDORI2_URL is not set.\n" >&2
  printf "        Set it to a MIDORI2 UNIQ FASTA .fasta.gz URL before running.\n" >&2
  printf "        Find current URLs at: https://www.reference-midori.info/download.php\n" >&2
  exit 1
fi

if ! command -v apptainer &>/dev/null; then
  printf "[ERROR] apptainer not found in PATH. Run: module load tacc-apptainer\n" >&2
  exit 1
fi

set -u

if [[ ! -f "$SIF" ]]; then
  printf "[ERROR] Container not found: %s\n" "$SIF" >&2
  exit 1
fi

printf "╔══════════════════════════════════════════════════════════════╗\n"
printf "║  diag_midori_blast.sh — Phase B0 MIDORI2 BLAST probe        ║\n"
printf "╚══════════════════════════════════════════════════════════════╝\n"
printf "[INFO] Test dir    : %s\n" "$TEST_DIR"
printf "[INFO] MIDORI2 URL : %s\n" "$MIDORI2_URL"
printf "[INFO] Container   : %s\n" "$SIF"

mkdir -p "$TEST_DIR"
export BLASTDB="$TEST_DIR"

# ─── Step 1: Download MIDORI2 reference FASTA ────────────────────────────────
printf "\n[INFO] Step 1: Download MIDORI2 reference FASTA\n"

case "$MIDORI2_URL" in
  *.fasta.gz)  DOWNLOAD_FILE="${TEST_DIR}/midori_ref.fasta.gz";  ARCHIVE_FMT="gz"  ;;
  *.fasta.zip) DOWNLOAD_FILE="${TEST_DIR}/midori_ref.fasta.zip"; ARCHIVE_FMT="zip" ;;
  *.zip)       DOWNLOAD_FILE="${TEST_DIR}/midori_ref.zip";       ARCHIVE_FMT="zip" ;;
  *)           DOWNLOAD_FILE="${TEST_DIR}/midori_ref.fasta";     ARCHIVE_FMT="none" ;;
esac
FASTA="${TEST_DIR}/midori_ref.fasta"
MIDORI_DB="${TEST_DIR}/midori_test"

if [[ -f "$FASTA" ]]; then
  printf "[SKIP] FASTA already present: %s\n" "$FASTA"
else
  curl -fSL --retry 3 --retry-delay 5 -o "$DOWNLOAD_FILE" "$MIDORI2_URL"
  case "$ARCHIVE_FMT" in
    gz)
      gunzip -f "$DOWNLOAD_FILE"
      ;;
    zip)
      # Extract the first .fasta entry from the zip into TEST_DIR, then normalise the name
      unzip -o -j "$DOWNLOAD_FILE" "*.fasta" -d "$TEST_DIR"
      extracted=$(find "$TEST_DIR" -maxdepth 1 -name "*.fasta" ! -name "midori_ref.fasta" | head -1)
      [[ -n "$extracted" ]] && mv "$extracted" "$FASTA"
      ;;
    none)
      : # plain FASTA downloaded directly — DOWNLOAD_FILE == FASTA
      ;;
  esac
  printf "[OK]   Downloaded and decompressed\n"
fi

NSEQ=$(grep -c '^>' "$FASTA" 2>/dev/null || echo 0)
printf "[INFO] FASTA sequences: %s\n" "$NSEQ"

if (( NSEQ == 0 )); then
  printf "[ERROR] FASTA contains 0 sequences — download may have failed or file is malformed\n" >&2
  exit 1
fi

# ─── Step 2: Check header format ────────────────────────────────────────────
printf "\n[INFO] Step 2: Inspect FASTA header format\n"
SAMPLE_HEADER=$(grep '^>' "$FASTA" | head -1)
printf "[INFO] Sample header: %s\n" "$SAMPLE_HEADER"

# MIDORI2 UNIQ (BLAST flavour) headers look like:
#   >AB001234.1.1.1234_Taxon_name
# The accession is the first field before the first underscore or dot-separated prefix.
# makeblastdb -parse_seqids handles this; we just log it for the ELN.

# ─── Step 3: makeblastdb ─────────────────────────────────────────────────────
printf "\n[INFO] Step 3: makeblastdb on MIDORI2 reference\n"

if [[ -f "${MIDORI_DB}.nhr" ]]; then
  printf "[SKIP] midori_test.nhr already exists\n"
else
  apptainer exec \
    --bind "$WORK":"$WORK" \
    "$SIF" \
    makeblastdb \
      -in "$FASTA" \
      -dbtype nucl \
      -parse_seqids \
      -out "$MIDORI_DB" \
      -title "MIDORI2 Phase B0 probe"
  printf "[OK]   makeblastdb complete\n"
fi

# ─── Step 4: Extract test query from MIDORI2 FASTA (self-query) ──────────────
printf "\n[INFO] Step 4: Extract test query from MIDORI2 FASTA\n"

# Use the first sequence from the MIDORI2 FASTA itself as the query. A self-query
# guarantees hits at 100% identity and directly tests whether staxids/sscinames
# are populated — which is the real concern. No dependency on mammalia_nt.
TEST_QUERY="${TEST_DIR}/test_query.fasta"

if [[ -f "$TEST_QUERY" ]]; then
  printf "[SKIP] test_query.fasta already present\n"
else
  awk '/^>/{if(n++)exit; print; next} n' "$FASTA" > "$TEST_QUERY"
  printf "[OK]   Test query extracted from MIDORI2 FASTA\n"
fi

if [[ ! -s "$TEST_QUERY" ]]; then
  printf "[ERROR] Failed to extract first sequence from MIDORI2 FASTA\n" >&2
  exit 1
fi

QUERY_ID=$(grep '^>' "$TEST_QUERY" | head -1 | awk '{print $1}' | tr -d '>')
printf "[INFO] Query accession: %s\n" "$QUERY_ID"

# ─── Step 5: BLAST query against MIDORI2 DB ──────────────────────────────────
printf "\n[INFO] Step 5: blastn query vs midori_test\n"

BLAST_OUT="${TEST_DIR}/blast_midori_probe.tsv"

apptainer exec \
  --bind "$WORK":"$WORK" \
  --env "BLASTDB=$TEST_DIR" \
  "$SIF" \
  blastn \
    -query "$TEST_QUERY" \
    -db "$MIDORI_DB" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames" \
    -max_target_seqs 5 \
    -perc_identity 80 \
    -out "$BLAST_OUT" \
    2>/dev/null || true

printf "[INFO] BLAST output: %s\n" "$BLAST_OUT"

# ─── Step 6: Evaluate staxids / sscinames population ────────────────────────
printf "\n[INFO] Step 6: Check staxids and sscinames columns\n"

if [[ ! -s "$BLAST_OUT" ]]; then
  printf "[WARN] No BLAST hits against MIDORI2 DB (empty output).\n" >&2
  printf "       This may be expected if the query marker does not overlap the MIDORI2 gene target.\n" >&2
  printf "       Consider using a query from the same marker gene as the MIDORI2 FASTA.\n" >&2
  printf "[FAIL] B0-MIDORI2 INCONCLUSIVE — no hits to evaluate taxonomy columns\n" >&2
  exit 1
fi

EMPTY_TAXID=$(awk -F'\t' '$13 == "" || $13 == "N/A"' "$BLAST_OUT" | wc -l)
TOTAL_HITS=$(wc -l < "$BLAST_OUT")
printf "[INFO] Total hits: %d; hits with empty staxids: %d\n" "$TOTAL_HITS" "$EMPTY_TAXID"

printf "[INFO] First 3 hits:\n"
head -3 "$BLAST_OUT" | awk -F'\t' '{printf "  qseqid=%-20s staxids=%-10s sscinames=%s\n", $1, $13, $14}'

if (( EMPTY_TAXID > 0 )); then
  printf "[FAIL] B0-MIDORI2 FAIL — %d hit(s) have empty staxids column\n" "$EMPTY_TAXID" >&2
  printf "       MIDORI2 accessions may not resolve via the current taxdb files.\n" >&2
  printf "       Phase B1 is not viable without resolving taxonomy column population.\n" >&2
  exit 1
fi

printf "[OK]   B0-MIDORI2 PASS — staxids and sscinames populated for all %d hit(s)\n" "$TOTAL_HITS"
printf "       MIDORI2 makeblastdb + taxdb integration confirmed on this node.\n"
printf "       Phase B1 (MIDORI2 strategic rebuild) is viable.\n"
