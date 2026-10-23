#!/usr/bin/env bash
set -eo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# scripts/smoke_test.sh  — end-to-end PLUMBING gate for the eDNA pipeline
#
# Purpose: prove that the pipeline runs from end to end and that every phase
#          boundary emits a structurally valid, non-empty artifact — then that
#          Phase 4 produced at least one real species-level call.
#
# This sits one rung BELOW scripts/validate_pipeline.sh. That script measures
# biological accuracy against a NanoSim mock community and needs simulated reads
# plus a ground-truth table. This one needs neither and answers the cheaper,
# more fundamental question: does the plumbing work at all?
#
# Capability / limitation (paired, per project convention):
#   DOES verify that each phase produces the artifact the next phase consumes,
#   that the Phase 4 column arithmetic is intact, and that taxonomy actually
#   resolves. Every blocker found in the 2026-08-04 review would have been caught
#   by the PRE-FLIGHT block below, which is where most of this script's value is.
#   Does NOT measure biological accuracy, sensitivity, or specificity — there is
#   no ground truth here. Passing this does NOT license reporting the taxonomy.
#
# Usage:
#   bash scripts/smoke_test.sh <raw_dir> <work_base> <db_dir> \
#        [threads] [recluster] [primer_file] [mode]
#
#   <raw_dir>    raw FASTQ tree (ignored when MODULES excludes Phase1)
#   <work_base>  scratch base; trimmed/ and consensus/ live under it
#   <db_dir>     directory holding the BLAST DB (stem comes from blast.db_name)
#
# Env:
#   SIF              container image      [default: Apptainer/edna-pipeline.sif]
#   PIPELINE_CONFIG  config path          [default: config/pipeline.yaml]
#   MODULES          phases to run/assert [default: all]
#   SKIP_RUN=1       assert against an EXISTING work_base without re-running.
#                    Lets the same gate grade a finished TACC batch run post-hoc.
#   BIND_PATH        extra path to bind into the container [default: work_base]
#
# Exit: 0 all assertions passed; 1 otherwise. Pre-flight failures abort
#       immediately (before burning compute); post-run assertions are collected
#       so one run reports every problem rather than only the first.
# ─────────────────────────────────────────────────────────────────────────────

RAW_DIR="${1:?Usage: smoke_test.sh <raw_dir> <work_base> <db_dir> [threads] [recluster] [primer_file] [mode]}"
WORK_BASE="${2:?need <work_base>}"
DB_DIR="${3:?need <db_dir>}"
THREADS="${4:-4}"
RECLUSTER="${5:-false}"
PRIMER_FILE="${6:-}"
MODE="${7:-otu}"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPTS_DIR/.." && pwd)"
SIF="${SIF:-$REPO_ROOT/Apptainer/edna-pipeline.sif}"
PIPELINE_CONFIG="${PIPELINE_CONFIG:-$REPO_ROOT/config/pipeline.yaml}"
MODULES="${MODULES:-all}"
BIND_PATH="${BIND_PATH:-$WORK_BASE}"

TRIMMED_DIR="$WORK_BASE/trimmed"
CONS_BASE="$WORK_BASE/consensus"
TAX_DIR="$CONS_BASE/combined/taxonomy"

FAILURES=0
ok()   { printf "[OK]   %s\n" "$*"; }
warn() { printf "[WARN] %s\n" "$*" >&2; }
fail() { printf "[FAIL] %s\n" "$*" >&2; FAILURES=$((FAILURES+1)); }
die()  { printf "[ERROR] %s\n" "$*" >&2; exit 1; }

wants() { [[ "$MODULES" == *all* || "$MODULES" == *"$1"* ]]; }

# non-empty file assertion used throughout
need_file() {
  local path="$1" what="$2"
  if [[ -s "$path" ]]; then ok "$what: $path"; return 0; fi
  fail "$what missing or empty: $path"; return 1
}

printf "╔══════════════════════════════════════════════════════════════╗\n"
printf "║  smoke_test.sh — pipeline plumbing gate                       ║\n"
printf "╚══════════════════════════════════════════════════════════════╝\n"
printf "[INFO] work_base=%s db_dir=%s modules=%s mode=%s recluster=%s\n" \
  "$WORK_BASE" "$DB_DIR" "$MODULES" "$MODE" "$RECLUSTER"
printf "[INFO] config=%s\n[INFO] sif=%s\n\n" "$PIPELINE_CONFIG" "$SIF"

# ═════════════════════════════════════════════════════════════════════════════
# PRE-FLIGHT — fail fast, before any compute is spent
# ═════════════════════════════════════════════════════════════════════════════
printf "── Pre-flight ──────────────────────────────────────────────────\n"

command -v apptainer &>/dev/null || die "apptainer not found (module load tacc-apptainer)"
[[ -f "$SIF" ]] || die "container image not found: $SIF"
[[ -f "$PIPELINE_CONFIG" ]] || die "config not found: $PIPELINE_CONFIG"
ok "apptainer + container image + config present"

# --- config keys resolve (a dropped key becomes an unbound-variable death mid-run)
# shellcheck source=lib/parse_yaml.sh
source "$SCRIPTS_DIR/lib/parse_yaml.sh"
eval "$(parse_yaml "$PIPELINE_CONFIG" CFG_)"
for k in CFG_blast_db_name CFG_blast_perc_identity CFG_blast_qcov_hsp_perc \
         CFG_blast_max_target_seqs CFG_blast_evalue CFG_blast_taxid \
         CFG_cutadapt_min_qual CFG_cutadapt_min_len CFG_cutadapt_max_len; do
  [[ -n "${!k:-}" ]] || die "config key empty or missing: $k (in $PIPELINE_CONFIG)"
done
ok "config keys resolve (db_name=$CFG_blast_db_name, pident=$CFG_blast_perc_identity, qcov=$CFG_blast_qcov_hsp_perc)"

# --- CRLF: a \r in any of these is fatal under Linux, and silently so
for f in "$SCRIPTS_DIR/pipeline.sh" "$SCRIPTS_DIR/lib/parse_yaml.sh" \
         "$SCRIPTS_DIR/lca_assignment.R" "$SCRIPTS_DIR/consensus.pl" \
         ${PRIMER_FILE:+"$PRIMER_FILE"}; do
  [[ -f "$f" ]] || continue
  if grep -q $'\r' "$f"; then
    die "CRLF line endings in $f — fatal under Linux. Fix: git add --renormalize . (see .gitattributes)"
  fi
done
ok "no CRLF in executed scripts or primer file"

# --- container tools. Rscript and ktImportText are the ones that historically
#     were absent from the shipped .sif while the .def claimed otherwise.
REQUIRED_TOOLS=(cutadapt fastqc NanoPlot multiqc blastn makeblastdb NGSpeciesID
                mafft cd-hit-est seqtk perl Rscript ktImportText)
[[ "$MODE" == "asv" ]] && REQUIRED_TOOLS+=(vsearch)
MISSING=""
for t in "${REQUIRED_TOOLS[@]}"; do
  apptainer exec "$SIF" bash -c "command -v '$t'" &>/dev/null || MISSING="$MISSING $t"
done
[[ -z "$MISSING" ]] || die "tools missing from container:$MISSING — rebuild from Apptainer/edna-pipeline.def"
ok "all ${#REQUIRED_TOOLS[@]} required tools present in container"

# --- R packages. tidyr is easy to miss: pipeline.sh guards the Rscript call with
#     `|| [WARN]`, so a missing package means taxonomy.tsv is never written and
#     the run still reports success.
apptainer exec "$SIF" Rscript -e 'suppressPackageStartupMessages({library(dplyr);library(tidyr)})' &>/dev/null \
  || die "R packages dplyr/tidyr not available in container — lca_assignment.R would fail SILENTLY"
ok "R packages dplyr + tidyr available"

# --- blastn argument compatibility.
#     blastn declares -taxids incompatible with -negative_taxids; passing both
#     aborts at argument-parse time on every database. This probe is the single
#     highest-value check in this script.
if grep -nE '^\s*FILTER_OPT=.*-taxids.*-negative_taxids' "$SCRIPTS_DIR/pipeline.sh" >/dev/null 2>&1; then
  die "pipeline.sh passes -taxids AND -negative_taxids to blastn — mutually incompatible, blastn will reject every call"
fi
#     Check the VERSION first: it gives a precise diagnosis where the argument
#     probe alone only says "rejected". `-taxids` arrived in BLAST+ 2.8, and a
#     bare `blast` conda spec under strict channel priority can silently resolve
#     to conda-forge's 2.5.0 (2018) — which fails every Phase 4 call.
BLASTN_VER="$(apptainer exec "$SIF" blastn -version 2>/dev/null | head -1 | sed 's/^blastn: *//; s/+.*//')"
BLASTN_MAJ="${BLASTN_VER%%.*}"
BLASTN_MIN="$(printf '%s' "$BLASTN_VER" | cut -d. -f2)"
if [[ -z "$BLASTN_VER" ]]; then
  die "could not determine blastn version in $SIF"
fi
if (( BLASTN_MAJ < 2 || (BLASTN_MAJ == 2 && BLASTN_MIN < 8) )); then
  die "blastn $BLASTN_VER is too old for -taxids (need >= 2.8).
       Rebuild the container: the .def must pin 'bioconda::blast>=2.15', because
       strict channel_priority otherwise resolves 'blast' to conda-forge's 2.5.0."
fi
ok "blastn $BLASTN_VER supports -taxids"

ARG_PROBE="$(apptainer exec "$SIF" blastn -task megablast \
              -db /nonexistent_db_probe -taxids "$CFG_blast_taxid" \
              -query /dev/null 2>&1 || true)"
if grep -q 'CArgException\|^USAGE' <<<"$ARG_PROBE"; then
  # Surface the actual Error lines, not the USAGE banner that precedes them.
  die "blastn rejected the configured filter arguments:
$(grep -E '^Error' <<<"$ARG_PROBE" | head -3 || tail -3 <<<"$ARG_PROBE")"
fi
ok "blastn accepts the configured filter arguments"

# --- database + taxonomy files
DB_PREFIX="$DB_DIR/$CFG_blast_db_name"
if compgen -G "${DB_PREFIX}*.n[hsi]r" >/dev/null || compgen -G "${DB_PREFIX}*.nal" >/dev/null \
   || compgen -G "${DB_PREFIX}*.nin" >/dev/null; then
  ok "BLAST DB present: $DB_PREFIX"
else
  die "no BLAST DB found at $DB_PREFIX.* — build one (scripts/build_smoke_db.sh or mammalia_nt_db.sh)"
fi
for tf in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do
  [[ -f "$DB_DIR/$tf" ]] || die "missing $DB_DIR/$tf — -taxids filtering and sscinames will not work"
done
ok "taxonomy lookup files present (taxdb.btd/.bti, taxonomy4blast.sqlite3)"

# --- DB format version. `-taxids` fails on a v4 database with
#     "Taxonomy filtering is not supported in v4 BLAST dbs". An old makeblastdb
#     (or one missing -blastdb_version 5) silently produces v4.
DB_INFO="$(apptainer exec --bind "$DB_DIR":"$DB_DIR" --env "BLASTDB=$DB_DIR" "$SIF" \
            blastdbcmd -db "$DB_PREFIX" -info 2>&1 || true)"
# blastdbcmd prints "BLASTDB Version: 5" (capital V) — match case-insensitively.
DB_VER="$(grep -oiE 'BLASTDB[[:space:]]+version:[[:space:]]*[0-9]+' <<<"$DB_INFO" | grep -oE '[0-9]+$' || true)"
if [[ -z "$DB_VER" ]]; then
  warn "could not read BLASTDB version from blastdbcmd -info (continuing)"
elif (( DB_VER < 5 )); then
  die "BLAST DB at $DB_PREFIX is version $DB_VER — -taxids requires v5.
       Rebuild it with 'makeblastdb -blastdb_version 5' (see scripts/build_smoke_db.sh)."
else
  ok "BLAST DB is version $DB_VER (v5 required for -taxids)"
fi

# --- refuse to run Phase 3 on top of a stale consensus tree. A directory full of
#     medaka_cl_id_* dirs with no consensus.fasta is an INCOMPLETE previous run;
#     the pipeline's own skip-guard does not catch it.
if [[ -z "${SKIP_RUN:-}" ]] && wants Phase3 && [[ -d "$CONS_BASE/per_barcode" ]]; then
  n_medaka=$(find "$CONS_BASE/per_barcode" -type d -name 'medaka_cl_id_*' 2>/dev/null | wc -l | tr -d ' ')
  n_cons=$(find "$CONS_BASE/per_barcode" -name 'consensus.fasta' 2>/dev/null | wc -l | tr -d ' ')
  if (( n_medaka > 0 && n_cons == 0 )); then
    die "CONS_BASE holds $n_medaka medaka_cl_id_* dirs but 0 consensus.fasta — an incomplete previous run. Use a virgin CONS_BASE."
  fi
fi
ok "pre-flight complete"

# ═════════════════════════════════════════════════════════════════════════════
# RUN
# ═════════════════════════════════════════════════════════════════════════════
if [[ -n "${SKIP_RUN:-}" ]]; then
  printf "\n[INFO] SKIP_RUN set — asserting against existing %s\n" "$WORK_BASE"
  WALL_STR="n/a (SKIP_RUN)"
else
  printf "\n── Running pipeline (modules=%s) ───────────────────────────────\n" "$MODULES"
  mkdir -p "$WORK_BASE"
  START=$(date +%s)
  apptainer exec \
    --bind "$REPO_ROOT":"$REPO_ROOT" --bind "$BIND_PATH":"$BIND_PATH" \
    --env "BLASTDB=$DB_DIR" --env "PIPELINE_CONFIG=$PIPELINE_CONFIG" \
    "$SIF" \
    bash "$SCRIPTS_DIR/pipeline.sh" \
      "$RAW_DIR" "$TRIMMED_DIR" "$CONS_BASE" "$DB_DIR" \
      "$THREADS" "$RECLUSTER" "$PRIMER_FILE" "$MODE" "$MODULES" \
    2>&1 | tee "$WORK_BASE/smoke_run.log"
  RC=${PIPESTATUS[0]}
  WALL=$(( $(date +%s) - START ))
  [[ $RC -eq 0 ]] || fail "pipeline exited $RC (see $WORK_BASE/smoke_run.log)"
  WALL_STR="${WALL}s"
  printf "[INFO] wall-time: %s\n" "$WALL_STR"
fi

# ═════════════════════════════════════════════════════════════════════════════
# ASSERTIONS — collected, so one run reports every problem
# ═════════════════════════════════════════════════════════════════════════════
printf "\n── Phase assertions ────────────────────────────────────────────\n"

if wants Phase1; then
  printf "[INFO] Phase 1 (trim)\n"
  mapfile -t TRIMMED < <(find "$TRIMMED_DIR" -name '*_trimmed.fastq.gz' 2>/dev/null | sort)
  if (( ${#TRIMMED[@]} == 0 )); then
    fail "no *_trimmed.fastq.gz under $TRIMMED_DIR"
  else
    ok "${#TRIMMED[@]} trimmed FASTQ produced"
    if [[ "$(zcat "${TRIMMED[0]}" 2>/dev/null | head -c1)" == "@" ]]; then
      ok "trimmed FASTQ starts with '@'"
    else
      fail "not valid FASTQ: ${TRIMMED[0]}"
    fi
  fi
fi

if wants Phase2; then
  printf "[INFO] Phase 2 (QC + merge)\n"
  mapfile -t MERGED < <(find "$TRIMMED_DIR/merged_fastq" -name '*_merged.fastq' 2>/dev/null | sort)
  if (( ${#MERGED[@]} == 0 )); then
    fail "no *_merged.fastq under $TRIMMED_DIR/merged_fastq"
  else
    ok "${#MERGED[@]} merged FASTQ produced"
    lines=$(wc -l < "${MERGED[0]}")
    (( lines % 4 == 0 )) && ok "merged FASTQ line count is a multiple of 4" \
                         || fail "truncated FASTQ (${lines} lines, not %4): ${MERGED[0]}"
  fi
fi

if wants Phase3; then
  printf "[INFO] Phase 3 (consensus)\n"
  n_cons=$(find "$CONS_BASE/per_barcode" -name 'consensus.fasta' -size +0 2>/dev/null | wc -l | tr -d ' ')
  # THE assertion whose absence let an incomplete run look successful.
  (( n_cons > 0 )) && ok "$n_cons non-empty per-barcode consensus.fasta" \
                   || fail "no non-empty medaka_cl_id_*/consensus.fasta — Phase 3 did not complete"
  if need_file "$CONS_BASE/combined/final_consensus.fasta" "final consensus"; then
    n=$(grep -c '^>' "$CONS_BASE/combined/final_consensus.fasta" || true)
    (( ${n:-0} > 0 )) && ok "final_consensus.fasta holds $n sequences" \
                      || fail "final_consensus.fasta has no FASTA records"
  fi
fi

if wants Phase4; then
  printf "[INFO] Phase 4 (taxonomy)\n"
  BLAST_ALL="$TAX_DIR/blast_all.tsv"
  if need_file "$BLAST_ALL" "BLAST raw output"; then
    # Column arithmetic: blast.tsv/blast_results.tsv are pure awk reshapes of
    # blast_all.tsv with no validation of their own. Spaces in sscinames/stitle
    # silently shift every downstream column if the awk loses -F'\t'.
    bad=$(awk -F'\t' 'NF!=17' "$BLAST_ALL" | wc -l | tr -d ' ')
    (( bad == 0 )) && ok "blast_all.tsv: all rows have 17 fields" \
                   || fail "blast_all.tsv: $bad row(s) not 17 fields"

    empty_tax=$(awk -F'\t' '$13=="" || $13=="0" || $13=="N/A"' "$BLAST_ALL" | wc -l | tr -d ' ')
    (( empty_tax == 0 )) && ok "staxids populated on every hit" \
                         || fail "$empty_tax hit(s) with empty staxids — DB lacks taxid mapping"

    empty_sci=$(awk -F'\t' '$14=="" || $14=="N/A"' "$BLAST_ALL" | wc -l | tr -d ' ')
    (( empty_sci == 0 )) && ok "sscinames populated on every hit" \
                         || fail "$empty_sci hit(s) with empty sscinames — taxdb.b* not found at query time"
  fi

  for spec in "blast.tsv:15" "blast_results.tsv:10"; do
    f="$TAX_DIR/${spec%%:*}"; want="${spec##*:}"
    if need_file "$f" "derived table"; then
      bad=$(awk -F'\t' -v w="$want" 'NF!=w' "$f" | wc -l | tr -d ' ')
      (( bad == 0 )) && ok "$(basename "$f"): all rows have $want fields" \
                     || fail "$(basename "$f"): $bad row(s) not $want fields (awk -F'\\t' lost?)"
    fi
  done

  TAX_TSV="$TAX_DIR/taxonomy.tsv"
  if need_file "$TAX_TSV" "LCA taxonomy"; then
    # On zero hits lca_assignment.R writes a header-only TSV and exits 0, so an
    # existence check alone reports a green run containing nothing.
    rows=$(( $(wc -l < "$TAX_TSV") - 1 ))
    (( rows > 0 )) && ok "taxonomy.tsv holds $rows assignment(s)" \
                   || fail "taxonomy.tsv is header-only — no assignments were made"

    if (( rows > 0 )); then
      n_species=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="lca_rank")r=i;if($i=="lca_taxon")t=i}next}
                              r&&t&&$r=="species"&&$t!="unresolved"&&$t!=""{n++} END{print n+0}' "$TAX_TSV")
      if (( n_species > 0 )); then
        ok "$n_species species-level call(s) — THE GATE"
        printf "[INFO] species called:\n"
        awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="lca_rank")r=i;if($i=="lca_taxon")t=i}next}
                    r&&t&&$r=="species"{print "         " $t}' "$TAX_TSV" | sort -u | head -20
      else
        fail "no species-level LCA call — taxonomy resolved nothing"
      fi
    fi
  fi

  KRONA="$TAX_DIR/krona.html"
  if [[ -s "$KRONA" ]] && grep -qi '<html' "$KRONA"; then
    ok "Krona chart rendered"
  else
    warn "Krona chart missing or not HTML: $KRONA (non-fatal)"
  fi
fi

# ═════════════════════════════════════════════════════════════════════════════
printf "\n── Summary ─────────────────────────────────────────────────────\n"
printf "  wall-time      : %s\n" "$WALL_STR"
printf "  work_base      : %s\n" "$WORK_BASE"
printf "  database       : %s\n" "$DB_PREFIX"
if (( FAILURES == 0 )); then
  printf "\n[OK]   SMOKE TEST PASSED — plumbing is sound.\n"
  printf "       This is NOT a validation of biological accuracy. For that, run\n"
  printf "       scripts/validate_pipeline.sh against a NanoSim mock community.\n"
  exit 0
else
  printf "\n[ERROR] SMOKE TEST FAILED — %d assertion(s) failed (see [FAIL] above).\n" "$FAILURES" >&2
  exit 1
fi
