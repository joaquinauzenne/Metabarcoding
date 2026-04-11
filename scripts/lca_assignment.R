#!/usr/bin/env Rscript
# lca_assignment.R — Lowest Common Ancestor assignment for BLAST results
#
# Reads a 10-column BLAST tabular output (format 6 with columns:
#   qseqid sseqid pident qcovs evalue bitscore staxids sscinames scomnames stitle)
# and assigns taxonomy at the lowest resolved level (species, genus, or family_or_higher)
# using a top-hit LCA approach.
#
# Usage: Rscript lca_assignment.R blast_results.tsv output_taxonomy.tsv
#
# Thresholds (configurable at top of file):
#   MIN_PIDENT  — minimum percent identity (default: 97)
#   MIN_QCOV    — minimum query coverage  (default: 90)
#   BITSCORE_MARGIN — keep hits within this fraction of the best bitscore (default: 0.98)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript lca_assignment.R <blast_results.tsv> <output_taxonomy.tsv>",
       call. = FALSE)
}
blast_file  <- args[1]
output_file <- args[2]

# ── Configurable thresholds ───────────────────────────────────────────────────
MIN_PIDENT      <- 97    # minimum percent identity
MIN_QCOV        <- 90    # minimum query coverage
BITSCORE_MARGIN <- 0.98  # keep top hits within 2% of best bitscore

# ── Read BLAST output (10-column CLI format) ──────────────────────────────────
blast <- tryCatch(
  read.delim(blast_file, header = FALSE,
    col.names = c("qseqid", "sseqid", "pident", "qcovs", "evalue",
                  "bitscore", "staxids", "sscinames", "scomnames", "stitle"),
    stringsAsFactors = FALSE),
  error = function(e) stop("Cannot read BLAST file: ", blast_file, "\n", e$message)
)

if (nrow(blast) == 0) {
  message("No BLAST hits found in ", blast_file, " — writing empty output.")
  write.table(
    data.frame(qseqid = character(0), lca_taxon = character(0),
               lca_rank = character(0), confidence = character(0)),
    output_file, sep = "\t", row.names = FALSE, quote = FALSE
  )
  quit(status = 0)
}

# ── Filter by identity and coverage ───────────────────────────────────────────
blast_filtered <- blast |>
  filter(pident >= MIN_PIDENT, qcovs >= MIN_QCOV) |>
  # Handle semicolon-delimited multi-taxon fields (when BLAST hits multiple taxa)
  separate_rows(sscinames, sep = ";") |>
  mutate(sscinames = trimws(sscinames)) |>
  # Drop environmental/uncultured sequences
  filter(!grepl("uncultured|environmental|synthetic|unclassified", sscinames,
                ignore.case = TRUE))

n_input    <- n_distinct(blast$qseqid)
n_filtered <- n_distinct(blast_filtered$qseqid)
message(sprintf("Filtered: %d queries → %d with pident≥%g%% and qcovs≥%g%%",
                n_input, n_filtered, MIN_PIDENT, MIN_QCOV))

# ── LCA assignment ────────────────────────────────────────────────────────────
lca_results <- blast_filtered |>
  group_by(qseqid) |>
  # Keep only hits within BITSCORE_MARGIN of the best bitscore for this query
  filter(bitscore >= max(bitscore) * BITSCORE_MARGIN) |>
  summarize(
    n_top_hits   = n(),
    n_species    = n_distinct(sscinames),
    best_pident  = max(pident),
    best_qcov    = max(qcovs),
    best_bitscore = max(bitscore),
    best_evalue  = min(evalue),
    best_species = if (n_distinct(sscinames) == 1) first(sscinames) else NA_character_,
    all_species  = paste(unique(sscinames), collapse = "; "),
    # Extract genus (first word of binomial name)
    genera       = paste(unique(sub(" .*", "", sscinames)), collapse = "; "),
    n_genera     = n_distinct(sub(" .*", "", sscinames)),
    common_names = paste(unique(scomnames[scomnames != "" & scomnames != "N/A"]),
                         collapse = "; "),
    .groups = "drop"
  ) |>
  mutate(
    # LCA rank: assign at species if all top hits agree, genus if same genus, else unresolved
    lca_rank = case_when(
      n_species == 1 ~ "species",
      n_genera  == 1 ~ "genus",
      TRUE           ~ "family_or_higher"
    ),
    lca_taxon = case_when(
      lca_rank == "species" ~ best_species,
      lca_rank == "genus"   ~ sub(";.*", "", genera),
      TRUE                  ~ "unresolved"
    ),
    # Confidence: based on identity and rank resolution
    confidence = case_when(
      best_pident >= 99 & lca_rank == "species"                   ~ "high",
      best_pident >= 97 & lca_rank %in% c("species", "genus")     ~ "medium",
      TRUE                                                         ~ "low"
    )
  )

message(sprintf("LCA assignment complete: %d queries → %d assigned (%d species, %d genus, %d unresolved)",
                n_filtered,
                nrow(lca_results),
                sum(lca_results$lca_rank == "species"),
                sum(lca_results$lca_rank == "genus"),
                sum(lca_results$lca_rank == "family_or_higher")))

# ── Write output ──────────────────────────────────────────────────────────────
write.table(lca_results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
message("Output written to: ", output_file)
