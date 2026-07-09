#!/usr/bin/env bash
set -eo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# scripts/simulate/simulate.sh  — M1 step 3 of 4
#
# Purpose: fetch per-species 12S/16S reference sequences (EDirect) and simulate
#          ONT reads from each with NanoSim (genome mode), one run per species so
#          ground truth is exact by construction. Reads are pooled into a single
#          simulated "barcode" FASTQ that pipeline.sh can consume unchanged.
#
# Why genome-mode-per-species (not Meta-NanoSim metagenome mode): it gives exact,
# unambiguous per-species read counts and avoids metagenome-mode's genome-list /
# abundance-file formats. Ground truth = the reads-per-species from mock_community.R.
#
# Capability / limitation (paired): DOES produce ONT-like reads with known source
# species for end-to-end validation. Does NOT reproduce PCR/primer bias, chimeras,
# or contamination — a passing validation bounds bioinformatic sequencing error.
#
# Dependencies:
#   - EDirect (esearch/efetch) on PATH, on the HOST (run scripts/diag_edirect_tacc.sh
#     first, or: export PATH="$HOME/edirect:$PATH"). Network → idev session.
#   - NanoSim inside the container's 'nanosim' env (see Apptainer/edna-pipeline.def).
#
# Usage (in an idev session, after mock_community.R and characterize.sh):
#   module load tacc-apptainer
#   bash scripts/simulate/simulate.sh <mock_dir> <out_raw_dir> [threads]
#     <mock_dir>     dir containing fetch_list.tsv + ground_truth.tsv
#     <out_raw_dir>  where the simulated barcode FASTQ is written (feeds pipeline.sh)
#
# Env: SIF, NANOSIM_ENV, MODEL_PREFIX (else read from nanosim_models/MODEL_PREFIX.txt),
#      SEED [42].
# ─────────────────────────────────────────────────────────────────────────────

MOCK_DIR="${1:?Usage: simulate.sh <mock_dir> <out_raw_dir> [threads]}"
OUT_RAW_DIR="${2:?Usage: simulate.sh <mock_dir> <out_raw_dir> [threads]}"
THREADS="${3:-8}"
SIF="${SIF:-Apptainer/edna-pipeline.sif}"
NANOSIM_ENV="${NANOSIM_ENV:-/opt/conda/envs/nanosim}"
MODELS_DIR="${MODELS_DIR:-${WORK:-$PWD}/Metabarcoding/nanosim_models}"
MODEL_PREFIX="${MODEL_PREFIX:-}"
SEED="${SEED:-42}"

FETCH_LIST="$MOCK_DIR/fetch_list.tsv"
REFS_DIR="$MOCK_DIR/refs"
SIM_DIR="$MOCK_DIR/nanosim_runs"
BARCODE_DIR="$OUT_RAW_DIR/sim_mock"
POOLED="$BARCODE_DIR/sim_mock.fastq"
FETCH_REPORT="$MOCK_DIR/fetch_report.tsv"

# ─── Guards ──────────────────────────────────────────────────────────────────
HOSTNAME_VAL="$(hostname)"
if [[ "$HOSTNAME_VAL" == login* ]] && [[ -z "${SLURM_JOB_ID:-}" ]]; then
  printf "[ERROR] Refusing to run on login node '%s'. Use idev.\n" "$HOSTNAME_VAL" >&2; exit 1
fi
[[ -f "$FETCH_LIST" ]] || { printf "[ERROR] fetch_list.tsv not found in %s (run mock_community.R first)\n" "$MOCK_DIR" >&2; exit 1; }
command -v apptainer &>/dev/null || { printf "[ERROR] apptainer not found. module load tacc-apptainer\n" >&2; exit 1; }
[[ -f "$SIF" ]] || { printf "[ERROR] Container not found: %s\n" "$SIF" >&2; exit 1; }
if ! command -v efetch &>/dev/null || ! command -v esearch &>/dev/null; then
  printf "[ERROR] EDirect (esearch/efetch) not on PATH.\n" >&2
  printf "        Run scripts/diag_edirect_tacc.sh, or: export PATH=\"\$HOME/edirect:\$PATH\"\n" >&2
  exit 1
fi

if [[ -z "$MODEL_PREFIX" ]]; then
  if [[ -f "$MODELS_DIR/MODEL_PREFIX.txt" ]]; then
    MODEL_PREFIX="$(cat "$MODELS_DIR/MODEL_PREFIX.txt")"
  else
    printf "[ERROR] No MODEL_PREFIX and no %s/MODEL_PREFIX.txt (run characterize.sh first)\n" "$MODELS_DIR" >&2
    exit 1
  fi
fi

# Bind the widest sensible path into the container so refs/model/output resolve.
BIND_PATH="${WORK:-$PWD}"

# ─── Preflight: confirm NanoSim genome-mode flags (do not assume) ────────────
printf "[INFO] Preflight: simulator.py genome --help\n"
HELP="$(apptainer exec "$SIF" "$NANOSIM_ENV/bin/simulator.py" genome --help 2>&1 || true)"
for flag in "-rg" "-c" "-o" "-n" "--fastq"; do
  if ! grep -q -- "$flag" <<<"$HELP"; then
    printf "[WARN] Expected flag '%s' not found in simulator.py genome --help.\n" "$flag" >&2
    printf "       NanoSim CLI may differ in this version — review the invocation below before trusting output.\n" >&2
  fi
done

mkdir -p "$REFS_DIR" "$SIM_DIR" "$BARCODE_DIR"
: > "$POOLED"
printf "scientific_name\tn_seqs_fetched\tref_len\tn_reads\tstatus\n" > "$FETCH_REPORT"

# ─── Per-species: fetch reference, simulate, relabel, pool ───────────────────
# skip header line of fetch_list.tsv
tail -n +2 "$FETCH_LIST" | while IFS=$'\t' read -r SP NREADS QUERY; do
  [[ -z "$SP" ]] && continue
  tag="$(printf '%s' "$SP" | tr ' /' '__')"
  ref_raw="$REFS_DIR/${tag}_all.fasta"
  ref="$REFS_DIR/${tag}.fasta"
  out_prefix="$SIM_DIR/$tag"

  printf "\n[INFO] === %s (n_reads=%s) ===\n" "$SP" "$NREADS"

  # Fetch (host EDirect). Guard SIGPIPE/failure; validate non-empty.
  esearch -db nuccore -query "$QUERY" 2>/dev/null | efetch -format fasta > "$ref_raw" 2>/dev/null || true
  nseq="$(grep -c '^>' "$ref_raw" 2>/dev/null || echo 0)"
  if (( nseq == 0 )); then
    printf "[WARN] No reference fetched for %s — skipping (will be a known false negative)\n" "$SP" >&2
    printf "%s\t0\t0\t%s\tNO_REFERENCE\n" "$SP" "$NREADS" >> "$FETCH_REPORT"
    continue
  fi

  # Pick the single longest record as the representative reference.
  awk '
    /^>/ { if (h!="" && length(seq)>maxlen){maxlen=length(seq);bh=h;bs=seq} h=$0; seq=""; next }
    { seq=seq $0 }
    END { if (h!="" && length(seq)>maxlen){bh=h;bs=seq} if(bh!=""){print bh; print bs} }
  ' "$ref_raw" > "$ref"
  ref_len="$(awk '!/^>/{n+=length($0)} END{print n+0}' "$ref")"
  printf "[INFO]   fetched %s seq(s); representative ref = %s bp\n" "$nseq" "$ref_len"

  # Simulate ONT reads from this species' reference.
  # VERIFY flags against 'simulator.py genome --help' (preflight above).
  apptainer exec --bind "$BIND_PATH":"$BIND_PATH" "$SIF" \
    "$NANOSIM_ENV/bin/simulator.py" genome \
      -rg "$ref" \
      -c "$MODEL_PREFIX" \
      -o "$out_prefix" \
      -n "$NREADS" \
      -t "$THREADS" \
      --fastq \
      --seed "$SEED" \
    || { printf "[WARN] NanoSim failed for %s — skipping\n" "$SP" >&2; \
         printf "%s\t%s\t%s\t%s\tSIM_FAILED\n" "$SP" "$nseq" "$ref_len" "$NREADS" >> "$FETCH_REPORT"; continue; }

  # Pool all reads from this run, relabeling headers with the source species tag
  # so the simulated FASTQ carries ground truth (readable by downstream audits).
  n_emitted=0
  while IFS= read -r -d '' rf; do
    awk -v tag="$tag" '
      NR%4==1 { c++; print "@SIM|" tag "|" c; next }
      { print }
    ' "$rf" >> "$POOLED"
    n_emitted=$(( n_emitted + $(( $(wc -l < "$rf") / 4 )) ))
  done < <(find "$SIM_DIR" -maxdepth 1 -name "${tag}*reads.fastq" -print0)

  printf "%s\t%s\t%s\t%s\tOK(%s reads)\n" "$SP" "$nseq" "$ref_len" "$NREADS" "$n_emitted" >> "$FETCH_REPORT"
done

if [[ ! -s "$POOLED" ]]; then
  printf "[ERROR] No reads simulated — see %s\n" "$FETCH_REPORT" >&2
  exit 1
fi

gzip -f "$POOLED"
printf "\n[OK]   Simulated barcode: %s.gz\n" "$POOLED"
printf "       Fetch/sim report:  %s\n" "$FETCH_REPORT"
printf "       Feed %s to pipeline.sh via scripts/validate_pipeline.sh\n" "$OUT_RAW_DIR"
