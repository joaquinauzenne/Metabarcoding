#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# scripts/simulate/mock_community.R  — M1 step 2 of 4
#
# Purpose: turn a USER-PROVIDED mammal checklist into a mock-community design for
#          NanoSim. Emits the ground-truth table and a per-species NCBI fetch list
#          that simulate.sh consumes. This script does NOT invent species and does
#          NOT touch the network — you supply the checklist; it only composes.
#
# Input checklist (CSV) — required column `scientific_name`; optional
#   `rel_abundance` (any positive numbers; normalized to sum to 1). Example:
#       scientific_name,rel_abundance
#       Alouatta seniculus,0.5
#       Dasypus novemcinctus,0.3
#       Cerdocyon thous,0.2
#   (Seed it from your Tatacoa/Magdalena mammal records. Mammal-only — no birds.)
#
# Usage:
#   Rscript scripts/simulate/mock_community.R <checklist.csv> <out_dir> [TOTAL_READS]
#
# Outputs (in <out_dir>):
#   ground_truth.tsv  — scientific_name, rel_abundance, n_reads  (the truth set)
#   fetch_list.tsv    — scientific_name, n_reads, esearch_query  (for simulate.sh)
# ─────────────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mock_community.R <checklist.csv> <out_dir> [TOTAL_READS]",
       call. = FALSE)
}
checklist_path <- args[1]
out_dir        <- args[2]
total_reads    <- if (length(args) >= 3) as.integer(args[3]) else 20000L

if (!file.exists(checklist_path)) {
  stop("Checklist not found: ", checklist_path, call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cl <- utils::read.csv(checklist_path, stringsAsFactors = FALSE, check.names = FALSE)
names(cl) <- tolower(trimws(names(cl)))
if (!"scientific_name" %in% names(cl)) {
  stop("Checklist must have a 'scientific_name' column. Found: ",
       paste(names(cl), collapse = ", "), call. = FALSE)
}

cl$scientific_name <- trimws(cl$scientific_name)
cl <- cl[nzchar(cl$scientific_name), , drop = FALSE]
cl <- cl[!duplicated(cl$scientific_name), , drop = FALSE]
if (nrow(cl) == 0) stop("Checklist has no usable species rows.", call. = FALSE)

# Relative abundances: use provided column, else even split.
if ("rel_abundance" %in% names(cl)) {
  ab <- suppressWarnings(as.numeric(cl$rel_abundance))
  if (anyNA(ab) || any(ab <= 0)) {
    stop("rel_abundance must be positive numbers for every row.", call. = FALSE)
  }
} else {
  ab <- rep(1, nrow(cl))
  message("No rel_abundance column — using an even split across species.")
}
rel <- ab / sum(ab)
n_reads <- as.integer(round(rel * total_reads))
# guarantee >=1 read per species so every truth species is actually simulated
n_reads[n_reads < 1L] <- 1L

# Per-species NCBI query: mitochondrial 12S/16S, sensible length window.
# simulate.sh runs this via EDirect (esearch | efetch). Adjust the markers/length
# window here if you target a different fragment.
esearch_query <- sprintf(
  '"%s"[Organism] AND (12S[Title] OR 16S[Title]) AND mitochondrion[Filter] AND 100:5000[SLEN]',
  cl$scientific_name
)

ground_truth <- data.frame(
  scientific_name = cl$scientific_name,
  rel_abundance   = round(rel, 6),
  n_reads         = n_reads,
  stringsAsFactors = FALSE
)
fetch_list <- data.frame(
  scientific_name = cl$scientific_name,
  n_reads         = n_reads,
  esearch_query   = esearch_query,
  stringsAsFactors = FALSE
)

utils::write.table(ground_truth, file.path(out_dir, "ground_truth.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(fetch_list, file.path(out_dir, "fetch_list.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("Mock community: %d species, %d total reads -> %s",
                nrow(ground_truth), sum(n_reads), out_dir))
message("  ground_truth.tsv and fetch_list.tsv written. Next: scripts/simulate/simulate.sh")
