#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# scripts/db_coverage_audit.R  — M2 (reviewer M2b): reference coverage audit
#
# Purpose: for a user mammal checklist, count how many species have >=1 12S and
#          >=1 16S record in the MIDORI2 reference FASTAs. Produces the methods
#          table a reviewer will ask for ("what fraction of the regional fauna is
#          actually represented at these markers?").
#
# Method (heuristic, format-agnostic): a species is "present" in a marker file if
# any header in that file contains its binomial, matched tolerantly with either a
# SPACE or an UNDERSCORE separator, case-insensitively. This does NOT assume
# MIDORI2's exact header layout. Failure mode is a graceful UNDERCOUNT (a genuine
# record with an unexpected spelling is missed) — never a false "present". State
# this limitation when reporting the numbers.
#
# Usage:
#   Rscript scripts/db_coverage_audit.R <checklist.csv> <midori_12S.fasta> <midori_16S.fasta> <out.tsv>
#     checklist.csv    must have a `scientific_name` column (same as mock_community.R)
#     midori_12S.fasta / midori_16S.fasta  the per-marker FASTAs downloaded by
#                                          build_midori_db.sh (in blastdb/midori_build/)
#
# Output: <out.tsv> (per-species has_12S/has_16S/has_both) + a printed summary.
# ─────────────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript db_coverage_audit.R <checklist.csv> <12S.fasta> <16S.fasta> <out.tsv>",
       call. = FALSE)
}
checklist_path <- args[1]; f12s <- args[2]; f16s <- args[3]; out_path <- args[4]

for (p in c(checklist_path, f12s, f16s)) {
  if (!file.exists(p)) stop("Input not found: ", p, call. = FALSE)
}

cl <- utils::read.csv(checklist_path, stringsAsFactors = FALSE, check.names = FALSE)
names(cl) <- tolower(trimws(names(cl)))
if (!"scientific_name" %in% names(cl)) {
  stop("Checklist must have a 'scientific_name' column.", call. = FALSE)
}
species <- unique(trimws(cl$scientific_name))
species <- species[nzchar(species)]
if (!length(species)) stop("Checklist has no species.", call. = FALSE)

# Read only header lines (memory-light) and lowercase once.
read_headers <- function(path) {
  con <- file(path, "r"); on.exit(close(con))
  hdrs <- character(0)
  repeat {
    chunk <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(chunk)) break
    h <- chunk[startsWith(chunk, ">")]
    if (length(h)) hdrs <- c(hdrs, h)
  }
  tolower(hdrs)
}

message("Reading headers...")
h12 <- read_headers(f12s)
h16 <- read_headers(f16s)
message(sprintf("  12S headers: %d ; 16S headers: %d", length(h12), length(h16)))

present_in <- function(sci, headers) {
  s <- tolower(trimws(sci))
  pat_space <- s
  pat_under <- gsub(" ", "_", s, fixed = TRUE)
  any(grepl(pat_space, headers, fixed = TRUE)) ||
    (pat_under != pat_space && any(grepl(pat_under, headers, fixed = TRUE)))
}

has_12s <- vapply(species, present_in, logical(1), headers = h12)
has_16s <- vapply(species, present_in, logical(1), headers = h16)

audit <- data.frame(
  scientific_name = species,
  has_12S = has_12s,
  has_16S = has_16s,
  has_both = has_12s & has_16s,
  has_either = has_12s | has_16s,
  stringsAsFactors = FALSE
)
utils::write.table(audit, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

n <- length(species)
pct <- function(x) sprintf("%d/%d (%.0f%%)", x, n, 100 * x / n)
cat("\n=== MIDORI2 reference coverage audit (heuristic; undercount-safe) ===\n")
cat(sprintf("checklist species:      %d\n", n))
cat(sprintf("with >=1 12S record:    %s\n", pct(sum(has_12s))))
cat(sprintf("with >=1 16S record:    %s\n", pct(sum(has_16s))))
cat(sprintf("with BOTH markers:      %s\n", pct(sum(audit$has_both))))
cat(sprintf("with NEITHER marker:    %s\n", pct(sum(!audit$has_either))))
missing <- species[!audit$has_either]
if (length(missing)) {
  cat(sprintf("\nNot found at either marker (%d): %s\n",
              length(missing), paste(missing, collapse = "; ")))
}
cat(sprintf("\nPer-species table: %s\n", out_path))
cat("NOTE: matching is a header-substring heuristic — treat counts as a lower bound.\n")
