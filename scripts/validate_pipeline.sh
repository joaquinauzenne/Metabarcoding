#!/usr/bin/env bash
set -eo pipefail
IFS=$'\n\t'

# ─────────────────────────────────────────────────────────────────────────────
# scripts/validate_pipeline.sh  — M1 step 4 of 4 (THE validation gate)
#
# Purpose: run the full pipeline (Phases 1-4) on the NanoSim-simulated reads, then
#          compare the combined LCA species calls to the mock-community ground
#          truth. Emits validation_report/ with sensitivity / precision / F1, the
#          false-negative / false-positive lists, and the run wall-time.
#
# This is the definition of "validated end-to-end" for the platform: a known
# mock community goes in, and we measure how faithfully species come back out —
# independent of real field data and (with a small ad-hoc DB) of the M2 build.
#
# Capability / limitation (paired): bounds the BIOINFORMATIC error (clustering,
# consensus, BLAST+LCA thresholds) against known truth. It does NOT validate
# wet-lab/PCR bias or real-sample contamination.
#
# Usage (idev + module load tacc-apptainer):
#   bash scripts/validate_pipeline.sh <sim_raw_dir> <work_base> <db_dir> <ground_truth.tsv> \
#        [threads] [mode] [primer_file]
#     <sim_raw_dir>       output of simulate.sh (contains sim_mock/*.fastq.gz)
#     <work_base>         scratch base for trimmed/consensus outputs
#     <db_dir>            BLAST DB dir (mammalia_midori once M2 is built, or a small ad-hoc DB)
#     <ground_truth.tsv>  from mock_community.R
# Env: SIF [Apptainer/edna-pipeline.sif]
# ─────────────────────────────────────────────────────────────────────────────

SIM_RAW_DIR="${1:?need <sim_raw_dir>}"
WORK_BASE="${2:?need <work_base>}"
DB_DIR="${3:?need <db_dir>}"
GROUND_TRUTH="${4:?need <ground_truth.tsv>}"
THREADS="${5:-8}"
MODE="${6:-otu}"
PRIMER_FILE="${7:-}"
SIF="${SIF:-Apptainer/edna-pipeline.sif}"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPTS_DIR/.." && pwd)"
TRIMMED_DIR="$WORK_BASE/trimmed"
CONS_BASE="$WORK_BASE/consensus"
REPORT="$WORK_BASE/validation_report"
TAX_TSV="$CONS_BASE/combined/taxonomy/taxonomy.tsv"
BIND_PATH="${WORK:-$REPO_ROOT}"

command -v apptainer &>/dev/null || { printf "[ERROR] apptainer not found. module load tacc-apptainer\n" >&2; exit 1; }
[[ -f "$SIF" ]] || { printf "[ERROR] Container not found: %s\n" "$SIF" >&2; exit 1; }
[[ -f "$GROUND_TRUTH" ]] || { printf "[ERROR] ground_truth.tsv not found: %s\n" "$GROUND_TRUTH" >&2; exit 1; }
mkdir -p "$REPORT"

# ─── Run the pipeline on simulated reads (timed) ─────────────────────────────
printf "[INFO] Running pipeline on simulated reads (mode=%s)\n" "$MODE"
START=$(date +%s)
apptainer exec --bind "$BIND_PATH":"$BIND_PATH" "$SIF" \
  bash "$SCRIPTS_DIR/pipeline.sh" \
    "$SIM_RAW_DIR" "$TRIMMED_DIR" "$CONS_BASE" "$DB_DIR" \
    "$THREADS" true "$PRIMER_FILE" "$MODE" all
END=$(date +%s)
WALL=$(( END - START ))
printf "%d seconds\n" "$WALL" > "$REPORT/walltime.txt"
printf "[OK]   Pipeline wall-time: %d s (recorded in %s/walltime.txt)\n" "$WALL" "$REPORT"

if [[ ! -s "$TAX_TSV" ]]; then
  printf "[ERROR] LCA taxonomy not produced: %s\n" "$TAX_TSV" >&2
  printf "        Phase 4 (or lca_assignment.R) did not complete — inspect the pipeline log.\n" >&2
  exit 1
fi

# ─── Compare LCA species calls to ground truth (base R) ──────────────────────
cat > "$REPORT/compare.R" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
truth_path <- args[1]; tax_path <- args[2]; out_dir <- args[3]

norm <- function(x) tolower(trimws(as.character(x)))

truth_df <- read.delim(truth_path, stringsAsFactors = FALSE)
truth <- unique(norm(truth_df$scientific_name))

tax <- read.delim(tax_path, stringsAsFactors = FALSE)
if (!all(c("lca_rank", "lca_taxon") %in% names(tax))) {
  stop("taxonomy.tsv missing lca_rank/lca_taxon columns; got: ",
       paste(names(tax), collapse = ", "))
}
called <- unique(norm(tax$lca_taxon[norm(tax$lca_rank) == "species"]))
called <- called[nzchar(called) & called != "unresolved"]

tp <- intersect(truth, called)
fn <- setdiff(truth, called)         # simulated but not recovered
fp <- setdiff(called, truth)         # called but not simulated

sens <- if (length(truth))  round(length(tp) / length(truth), 4)  else NA
prec <- if (length(called)) round(length(tp) / length(called), 4) else NA
f1   <- if (!is.na(sens) && !is.na(prec) && (sens + prec) > 0)
          round(2 * sens * prec / (sens + prec), 4) else NA

per <- data.frame(
  species   = sort(union(truth, called)),
  stringsAsFactors = FALSE
)
per$in_truth  <- per$species %in% truth
per$in_called <- per$species %in% called
per$class <- ifelse(per$in_truth & per$in_called, "TP",
              ifelse(per$in_truth & !per$in_called, "FN", "FP"))
write.table(per, file.path(out_dir, "per_species.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

lines <- c(
  "=== M1 validation: simulated mock community vs pipeline LCA calls ===",
  sprintf("truth species (simulated): %d", length(truth)),
  sprintf("called species (LCA):      %d", length(called)),
  sprintf("true positives:            %d", length(tp)),
  sprintf("sensitivity (recall):      %s", ifelse(is.na(sens), "NA", sens)),
  sprintf("precision:                 %s", ifelse(is.na(prec), "NA", prec)),
  sprintf("F1:                        %s", ifelse(is.na(f1), "NA", f1)),
  "",
  sprintf("false negatives (%d): %s", length(fn), paste(fn, collapse = "; ")),
  sprintf("false positives (%d): %s", length(fp), paste(fp, collapse = "; "))
)
writeLines(lines, file.path(out_dir, "metrics.txt"))
cat(paste(lines, collapse = "\n"), "\n")
RSCRIPT

printf "\n[INFO] Comparing LCA calls to ground truth\n"
apptainer exec --bind "$BIND_PATH":"$BIND_PATH" "$SIF" \
  Rscript "$REPORT/compare.R" "$GROUND_TRUTH" "$TAX_TSV" "$REPORT"

printf "\n[OK]   Validation report: %s\n" "$REPORT"
printf "       metrics.txt, per_species.tsv, walltime.txt\n"
printf "       Record the sensitivity/precision + wall-time in a new ELN entry.\n"
