#!/usr/bin/env bash
set -eo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# scripts/build_smoke_db.sh  — build the small local BLAST DB used by the
#                              end-to-end smoke test (scripts/smoke_test.sh)
#
# Purpose: turn an NCBI web-BLAST "Alignment Descriptions" CSV export into a
#          small, taxid-resolving BLAST database, so Phase 4 can be exercised
#          locally without the multi-GB MIDORI2 build and without TACC.
#
# Why seed from the web-BLAST export rather than a generic mammal panel: those
# accessions are the hits NCBI actually returned FOR THIS PROJECT'S OWN consensus
# sequences (they are how Phase 4 was done by hand before it was automated). A
# database built from them is near-guaranteed to produce hits, so an empty
# taxonomy.tsv means "the plumbing is broken", not "the reference missed".
#
# Capability / limitation (paired, per project convention):
#   DOES produce a real, taxid-resolving BLAST DB whose staxids/sscinames populate,
#   sufficient to prove Phase 4 end to end against known-good expectations.
#   Does NOT constitute a reference database for producing results. It is one
#   query's top hits — a few dozen species, heavily biased to Neotropical primates.
#   Sequences absent from it will correctly return nothing. Never publish findings
#   from this DB; use scripts/mammalia_nt_db.sh for that.
#
# NOT a substitute for mammalia_nt. The DB is deliberately named
# `smoke_neotrop` so it can never be mistaken for the real M2 reference build.
#
# Usage:
#   bash scripts/build_smoke_db.sh [alignment_descriptions.csv]
#     default CSV: examples/2ER9WHC1013-Alignment-Descriptions.csv
#
# Env:
#   OUT_DIR        DB output dir            [default: blastdb/smoke]
#   DB_NAME        makeblastdb -out stem    [default: smoke_neotrop]
#                  (MUST match blast.db_name in config/pipeline.smoke.yaml)
#   SIF            container image          [default: Apptainer/edna-pipeline.sif]
#   TAXDB_TARBALL  taxdb archive            [default: blastdb/taxdb.tar.gz]
#   NCBI_API_KEY   optional; raises the NCBI rate limit from 3/s to 10/s
#   VERIFY_TAXIDS  cross-check taxids against NCBI esummary [default: true]
# ─────────────────────────────────────────────────────────────────────────────

CSV="${1:-examples/2ER9WHC1013-Alignment-Descriptions.csv}"
OUT_DIR="${OUT_DIR:-blastdb/smoke}"
DB_NAME="${DB_NAME:-smoke_neotrop}"
SIF="${SIF:-Apptainer/edna-pipeline.sif}"
TAXDB_TARBALL="${TAXDB_TARBALL:-blastdb/taxdb.tar.gz}"
VERIFY_TAXIDS="${VERIFY_TAXIDS:-true}"

DB_PREFIX="$OUT_DIR/$DB_NAME"
BUILD_JSON="${DB_PREFIX}.build.json"
REF_FASTA="$OUT_DIR/ref.fasta"
TAXID_MAP="$OUT_DIR/acc_taxid.map"
ACC_LIST="$OUT_DIR/accessions.txt"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ─── Guards ──────────────────────────────────────────────────────────────────
[[ -f "$CSV" ]] || { printf "[ERROR] CSV not found: %s\n" "$CSV" >&2; exit 1; }
command -v curl &>/dev/null || { printf "[ERROR] curl not found\n" >&2; exit 1; }
command -v apptainer &>/dev/null || { printf "[ERROR] apptainer not found. module load tacc-apptainer\n" >&2; exit 1; }
[[ -f "$SIF" ]] || { printf "[ERROR] Container not found: %s\n" "$SIF" >&2; exit 1; }

PY=""
for c in python3 python; do command -v "$c" &>/dev/null && { PY="$c"; break; }; done
[[ -n "$PY" ]] || { printf "[ERROR] Need python3 (or python) to parse the CSV\n" >&2; exit 1; }

mkdir -p "$OUT_DIR"

# ─── Step 1: parse accession + taxid out of the web-BLAST CSV ────────────────
# The Accession column holds an Excel =HYPERLINK(...) formula (and its header is
# literally "Accession  ", with trailing spaces), so this is a real parse, not a
# cut -f. Fields may span lines, hence Python's csv module rather than awk.
printf "[INFO] === Step 1: parse accessions from %s ===\n" "$CSV"
"$PY" - "$CSV" "$TAXID_MAP" "$ACC_LIST" <<'PYEOF'
import csv, re, sys
csv_path, map_path, acc_path = sys.argv[1], sys.argv[2], sys.argv[3]
with open(csv_path, encoding="utf-8-sig", newline="") as fh:
    rows = list(csv.reader(fh))
if not rows:
    sys.exit("[ERROR] CSV is empty")
hdr = [h.strip() for h in rows[0]]
for col in ("Accession", "Taxid"):
    if col not in hdr:
        sys.exit(f"[ERROR] CSV missing '{col}' column; found: {hdr}")
ai, ti = hdr.index("Accession"), hdr.index("Taxid")
si = hdr.index("Scientific Name") if "Scientific Name" in hdr else None

pat = re.compile(r"([A-Z]{1,2}_?\d{5,}\.\d+)")
pairs, short, noacc, badtax = {}, 0, 0, 0
for r in rows[1:]:
    if len(r) <= max(ai, ti):
        short += 1; continue
    m = pat.findall(r[ai])
    if not m:
        noacc += 1; continue
    taxid = r[ti].strip()
    if not taxid.isdigit():
        badtax += 1; continue
    # last match: the formula contains the accession in both URL and label
    pairs[m[-1]] = (taxid, r[si].strip() if si is not None else "")

if not pairs:
    sys.exit("[ERROR] No usable accession/taxid pairs parsed")

with open(map_path, "w", newline="\n") as fh:
    for acc, (taxid, _) in sorted(pairs.items()):
        fh.write(f"{acc} {taxid}\n")
with open(acc_path, "w", newline="\n") as fh:
    for acc in sorted(pairs):
        fh.write(acc + "\n")

species = {v[1] for v in pairs.values() if v[1]}
print(f"[OK]   parsed {len(pairs)} accessions, {len({v[0] for v in pairs.values()})} taxids, {len(species)} species")
if short or noacc or badtax:
    print(f"[WARN] skipped rows -> short:{short} no-accession:{noacc} non-numeric-taxid:{badtax}")
PYEOF

N_ACC=$(wc -l < "$ACC_LIST" | tr -d ' ')
(( N_ACC > 0 )) || { printf "[ERROR] No accessions parsed\n" >&2; exit 1; }

# ─── Step 2: cross-check taxids against NCBI (catches a stale/edited CSV) ─────
if [[ "$VERIFY_TAXIDS" == "true" ]]; then
  printf "\n[INFO] === Step 2: verify taxids against NCBI esummary ===\n"
  VERIFY_OUT="$OUT_DIR/taxid_verify.tsv"; : > "$VERIFY_OUT"
  KEY_PARAM=""; [[ -n "${NCBI_API_KEY:-}" ]] && KEY_PARAM="&api_key=${NCBI_API_KEY}"
  batch=""; n=0; mismatch=0; checked=0
  flush_batch() {
    [[ -z "$batch" ]] && return 0
    local json
    json="$(curl -sS --retry 3 --retry-delay 2 \
      "${EUTILS}/esummary.fcgi?db=nuccore&id=${batch}&retmode=json${KEY_PARAM}" || true)"
    if [[ -n "$json" ]]; then
      # NOTE: no backslashes inside f-string expressions here — that is a syntax
      # error before Python 3.12, and it fails silently behind the `|| true`.
      printf '%s' "$json" | "$PY" -c '
import json, sys
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
res = d.get("result", {})
for uid in res.get("uids", []):
    r = res.get(uid, {})
    ver = r.get("accessionversion") or r.get("caption", "")
    tax = r.get("taxid", "")
    if ver and tax:
        print(str(ver) + "\t" + str(tax))
' >> "$VERIFY_OUT" || true
    fi
    batch=""
    sleep 0.4
  }
  while IFS= read -r acc; do
    batch="${batch:+$batch,}$acc"; n=$((n+1))
    (( n % 50 == 0 )) && flush_batch
  done < "$ACC_LIST"
  flush_batch

  if [[ -s "$VERIFY_OUT" ]]; then
    # compare NCBI's taxid to the CSV's for accessions we could look up
    while IFS=$'\t' read -r vacc vtax; do
      [[ -z "$vacc" || -z "$vtax" ]] && continue
      ctax="$(awk -v a="$vacc" '$1==a{print $2; exit}' "$TAXID_MAP")"
      [[ -z "$ctax" ]] && continue
      checked=$((checked+1))
      if [[ "$ctax" != "$vtax" ]]; then
        mismatch=$((mismatch+1))
        printf "[WARN] taxid mismatch %s: CSV=%s NCBI=%s (using NCBI)\n" "$vacc" "$ctax" "$vtax" >&2
        # NCBI is authoritative
        awk -v a="$vacc" -v t="$vtax" 'BEGIN{OFS=" "} $1==a{$2=t} {print}' "$TAXID_MAP" > "$TAXID_MAP.tmp" \
          && mv "$TAXID_MAP.tmp" "$TAXID_MAP"
      fi
    done < "$VERIFY_OUT"
    printf "[OK]   taxid cross-check: %d verified, %d corrected\n" "$checked" "$mismatch"
  else
    printf "[WARN] esummary returned nothing — skipping cross-check (offline?)\n" >&2
  fi
else
  printf "\n[INFO] === Step 2: taxid cross-check SKIPPED (VERIFY_TAXIDS=false) ===\n"
fi

# ─── Step 3: fetch the reference sequences ───────────────────────────────────
printf "\n[INFO] === Step 3: fetch %s sequences from NCBI ===\n" "$N_ACC"
if [[ -s "$REF_FASTA" ]]; then
  printf "[SKIP] %s already present\n" "$REF_FASTA"
else
  : > "$REF_FASTA"
  KEY_PARAM=""; [[ -n "${NCBI_API_KEY:-}" ]] && KEY_PARAM="&api_key=${NCBI_API_KEY}"
  batch=""; n=0
  fetch_batch() {
    [[ -z "$batch" ]] && return 0
    curl -sS --retry 3 --retry-delay 2 \
      "${EUTILS}/efetch.fcgi?db=nuccore&id=${batch}&rettype=fasta&retmode=text${KEY_PARAM}" \
      >> "$REF_FASTA" || printf "[WARN] efetch batch failed\n" >&2
    batch=""
    sleep 0.4
  }
  while IFS= read -r acc; do
    batch="${batch:+$batch,}$acc"; n=$((n+1))
    (( n % 50 == 0 )) && { fetch_batch; printf "[INFO]   fetched %d/%d\n" "$n" "$N_ACC"; }
  done < "$ACC_LIST"
  fetch_batch
fi

N_SEQ=$(grep -c '^>' "$REF_FASTA" 2>/dev/null || true)
N_SEQ="${N_SEQ:-0}"
printf "[INFO] sequences in %s: %s\n" "$REF_FASTA" "$N_SEQ"
(( N_SEQ > 0 )) || { printf "[ERROR] No sequences fetched — check network\n" >&2; exit 1; }
if (( N_SEQ < N_ACC )); then
  printf "[WARN] fetched %s of %s accessions; the DB will be correspondingly thinner\n" "$N_SEQ" "$N_ACC" >&2
fi

# ─── Step 4: taxonomy lookup files ───────────────────────────────────────────
# -taxids filtering and the staxids/sscinames columns require these next to the DB.
printf "\n[INFO] === Step 4: taxonomy files ===\n"
if [[ -f "$OUT_DIR/taxdb.btd" && -f "$OUT_DIR/taxonomy4blast.sqlite3" ]]; then
  printf "[SKIP] taxdb files already present in %s\n" "$OUT_DIR"
elif [[ -f "$TAXDB_TARBALL" ]]; then
  tar -xzf "$TAXDB_TARBALL" -C "$OUT_DIR"
  printf "[OK]   extracted %s -> %s\n" "$TAXDB_TARBALL" "$OUT_DIR"
else
  printf "[WARN] %s not found; downloading taxdb from NCBI\n" "$TAXDB_TARBALL" >&2
  curl -fSL --retry 3 -o "$OUT_DIR/taxdb.tar.gz" \
    "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz" \
    && tar -xzf "$OUT_DIR/taxdb.tar.gz" -C "$OUT_DIR"
fi
for tf in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do
  [[ -f "$OUT_DIR/$tf" ]] && printf "[OK]   %s\n" "$tf" \
    || printf "[WARN] missing %s — -taxids and sscinames will not work\n" "$tf" >&2
done

# ─── Step 5: makeblastdb ─────────────────────────────────────────────────────
printf "\n[INFO] === Step 5: makeblastdb -> %s ===\n" "$DB_PREFIX"
BIND_ROOT="$(cd "$(dirname "$OUT_DIR")" && pwd)"
apptainer exec --bind "$BIND_ROOT":"$BIND_ROOT" "$SIF" \
  makeblastdb \
    -in "$REF_FASTA" \
    -dbtype nucl \
    -parse_seqids \
    -blastdb_version 5 \
    -taxid_map "$TAXID_MAP" \
    -out "$DB_PREFIX" \
    -title "SMOKE-TEST DB (${N_SEQ} seqs from ${CSV}) - NOT a reference database"
printf "[OK]   makeblastdb complete\n"

# ─── Step 6: self-query verification ─────────────────────────────────────────
# Same criterion as build_midori_db.sh: BLAST the first sequence back and confirm
# the taxonomy columns actually populate (proves taxid_map + taxdb together).
printf "\n[INFO] === Step 6: verify taxonomy resolution (self-query) ===\n"
TEST_QUERY="$OUT_DIR/test_query.fasta"
awk '/^>/{if(n++)exit; print; next} n' "$REF_FASTA" > "$TEST_QUERY"
BLAST_OUT="$OUT_DIR/verify.tsv"
apptainer exec --bind "$BIND_ROOT":"$BIND_ROOT" --env "BLASTDB=$(cd "$OUT_DIR" && pwd)" "$SIF" \
  blastn -query "$TEST_QUERY" -db "$DB_PREFIX" \
    -outfmt "6 qseqid sseqid pident staxids sscinames" \
    -max_target_seqs 3 -perc_identity 90 -out "$BLAST_OUT" 2>/dev/null || true

VERIFY_STATUS="FAIL"
if [[ -s "$BLAST_OUT" ]]; then
  EMPTY_TAX=$(awk -F'\t' '$4=="" || $4=="N/A" || $4=="0"' "$BLAST_OUT" | wc -l | tr -d ' ')
  EMPTY_SCI=$(awk -F'\t' '$5=="" || $5=="N/A"' "$BLAST_OUT" | wc -l | tr -d ' ')
  if (( EMPTY_TAX == 0 && EMPTY_SCI == 0 )); then
    VERIFY_STATUS="PASS"
    printf "[OK]   taxonomy columns populated. First hit: %s\n" "$(head -1 "$BLAST_OUT")"
  else
    printf "[WARN] %s hit(s) empty staxids, %s empty sscinames — check taxdb/taxid_map\n" \
      "$EMPTY_TAX" "$EMPTY_SCI" >&2
  fi
else
  printf "[WARN] self-query returned no hits — verification inconclusive\n" >&2
fi

# ─── Step 7: provenance ──────────────────────────────────────────────────────
sha() { sha256sum "$1" 2>/dev/null | awk '{print $1}'; }
BLAST_VER="$(apptainer exec "$SIF" blastn -version 2>/dev/null | head -1 || echo unknown)"
{
  printf '{\n'
  printf '  "db_name": "%s",\n' "$DB_NAME"
  printf '  "purpose": "smoke test only - NOT a reference database for results",\n'
  printf '  "source_csv": "%s",\n' "$CSV"
  printf '  "source_csv_sha256": "%s",\n' "$(sha "$CSV")"
  printf '  "built_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '  "blast_version": "%s",\n' "$BLAST_VER"
  printf '  "n_accessions_requested": %s,\n' "$N_ACC"
  printf '  "n_sequences_fetched": %s,\n' "$N_SEQ"
  printf '  "ref_fasta_sha256": "%s",\n' "$(sha "$REF_FASTA")"
  printf '  "taxid_map_sha256": "%s",\n' "$(sha "$TAXID_MAP")"
  printf '  "mammalia_scoping": "query-time -taxids 40674 (pipeline.sh); Hominidae dropped in lca_assignment.R",\n'
  printf '  "taxonomy_self_query": "%s"\n' "$VERIFY_STATUS"
  printf '}\n'
} > "$BUILD_JSON"

printf "\n[OK]   Smoke DB built: %s.*\n" "$DB_PREFIX"
printf "       Provenance:   %s\n" "$BUILD_JSON"
printf "       Verification: %s\n" "$VERIFY_STATUS"
printf "       Point the pipeline at it with DB_DIR=%s and blast.db_name=%s\n" "$OUT_DIR" "$DB_NAME"
printf "       (config/pipeline.smoke.yaml already sets that db_name.)\n"
[[ "$VERIFY_STATUS" == "PASS" ]] || {
  printf "[ERROR] Taxonomy verification did not PASS — Phase 4 would produce empty staxids.\n" >&2
  exit 1
}
