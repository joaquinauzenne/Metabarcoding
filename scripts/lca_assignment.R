#!/usr/bin/env Rscript
# lca_assignment.R — Lowest Common Ancestor assignment for BLAST results
#
# Reads a 10-column BLAST tabular output (format 6 with columns:
#   qseqid sseqid pident qcovs evalue bitscore staxids sscinames scomnames stitle)
# and assigns taxonomy at the lowest resolved level (species, genus, or family_or_higher)
# using a top-hit LCA approach.
#
# Usage: Rscript lca_assignment.R <blast_results.tsv> <output_taxonomy.tsv> [exclude_regex]
#
# Thresholds are read from config/pipeline.yaml (`blast:` block) so the config is
# the single source of truth; the literals below are only fallbacks for standalone
# runs where the config is unavailable. Override the config path with PIPELINE_CONFIG.
#   MIN_PIDENT      <- blast.perc_identity   (fallback: 97)
#   MIN_QCOV        <- blast.qcov_hsp_perc   (fallback: 90)
#   BITSCORE_MARGIN <- blast.bitscore_margin (fallback: 0.98)
#
# This script also applies the Hominidae exclusion that blastn cannot — see the
# "Exclude Hominidae" section below for why it lives here rather than in Phase 4.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript lca_assignment.R <blast_results.tsv> <output_taxonomy.tsv> [exclude_regex]",
       call. = FALSE)
}
blast_file  <- args[1]
output_file <- args[2]

# ── Thresholds: config-driven with literal fallbacks ──────────────────────────
MIN_PIDENT      <- 97    # minimum percent identity
MIN_QCOV        <- 90    # minimum query coverage
BITSCORE_MARGIN <- 0.98  # keep top hits within 2% of best bitscore

.script_dir <- {
  cargs <- commandArgs(trailingOnly = FALSE)
  farg  <- grep("^--file=", cargs, value = TRUE)
  if (length(farg) == 1L) dirname(normalizePath(sub("^--file=", "", farg))) else getwd()
}
config_path <- Sys.getenv("PIPELINE_CONFIG",
                          unset = file.path(.script_dir, "..", "config", "pipeline.yaml"))

if (file.exists(config_path) && requireNamespace("yaml", quietly = TRUE)) {
  .cfg <- yaml::read_yaml(config_path)$blast
  if (!is.null(.cfg$perc_identity))   MIN_PIDENT      <- as.numeric(.cfg$perc_identity)
  if (!is.null(.cfg$qcov_hsp_perc))   MIN_QCOV        <- as.numeric(.cfg$qcov_hsp_perc)
  if (!is.null(.cfg$bitscore_margin)) BITSCORE_MARGIN <- as.numeric(.cfg$bitscore_margin)
  message(sprintf("Thresholds from %s: pident>=%g, qcov>=%g, bitscore_margin=%g",
                  config_path, MIN_PIDENT, MIN_QCOV, BITSCORE_MARGIN))
} else {
  message(sprintf("Config not readable (%s) - using built-in defaults: pident>=%g, qcov>=%g, bitscore_margin=%g",
                  config_path, MIN_PIDENT, MIN_QCOV, BITSCORE_MARGIN))
}

# Hominidae genera. Overridable via 3rd arg or LCA_EXCLUDE_REGEX; set to "" to disable.
EXCLUDE_PATTERN <- if (length(args) >= 3L) args[3] else
  Sys.getenv("LCA_EXCLUDE_REGEX", unset = "\\b(Homo|Pan|Gorilla|Pongo)\\b")

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

# ── Exclude Hominidae (human / great-ape contamination) ───────────────────────
# This filter lives HERE, not in the blastn call, because blastn declares -taxids
# and -negative_taxids mutually incompatible: passing both aborts at argument-parse
# time before the database is opened. Phase 4 therefore scopes to Mammalia with a
# positive `-taxids 40674`, and the Hominidae exclusion is applied at this step.
#
# Matching is on the four Hominidae genera rather than on taxid 9604, because the
# staxids column alone cannot be resolved to descendants of 9604 without walking
# the NCBI taxonomy tree. Override with a 3rd arg or LCA_EXCLUDE_REGEX.
#
# SEMANTICS: this runs BEFORE the bitscore-margin filter below. A query whose single
# best hit is human but whose next-best hit is a genuine mammal will therefore
# resolve to that mammal rather than being discarded outright. That is intended.
#
# We COUNT before dropping: in eDNA, human read-through is a reportable
# contamination metric, not merely noise to discard silently.
if (nzchar(EXCLUDE_PATTERN)) {
  is_hominid <- grepl(EXCLUDE_PATTERN, blast_filtered$sscinames, ignore.case = TRUE)
  n_hom_hits    <- sum(is_hominid)
  n_hom_queries <- n_distinct(blast_filtered$qseqid[is_hominid])
  hom_taxa      <- sort(unique(blast_filtered$sscinames[is_hominid]))
  if (n_hom_hits > 0) {
    message(sprintf("Hominidae excluded: %d hit(s) across %d quer(ies) [%s]",
                    n_hom_hits, n_hom_queries, paste(hom_taxa, collapse = "; ")))
  } else {
    message("Hominidae excluded: 0 hits (no human/great-ape contamination detected)")
  }
  blast_filtered <- blast_filtered[!is_hominid, , drop = FALSE]
} else {
  message("Hominidae exclusion DISABLED (empty exclude pattern)")
}

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
