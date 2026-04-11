#!/usr/bin/env Rscript
# asv_otu_comparison.R — Compare taxonomic assignments across pipeline branches
#
# Reads LCA taxonomy TSVs from three pipeline branches and produces a
# species concordance matrix to identify discordant assignments:
#   - concordant: all three pipelines agree (high confidence)
#   - asv_only:   NGSpeciesID-ASV detected but OTU did not (ASV adds resolution)
#   - dada2_only: DADA2 detected but NGSpeciesID-ASV did not (verify manually)
#   - otu_only:   OTU detected but ASV did not (possible OTU artifact)
#
# Usage:
#   Rscript asv_otu_comparison.R <otu_tax.tsv> <asv_tax.tsv> <dada2_tax.tsv> <output_dir>
#
# Input TSVs must be produced by lca_assignment.R and contain columns:
#   qseqid, lca_taxon, lca_rank, confidence, all_species, ...
#
# Outputs (in output_dir/):
#   species_concordance.csv   — species × pipeline detection matrix with flags
#   concordance_summary.txt   — printed summary counts

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript asv_otu_comparison.R <otu_tax.tsv> <asv_tax.tsv> <dada2_tax.tsv> <output_dir>",
       call. = FALSE)
}

otu_file    <- args[1]
asv_file    <- args[2]
dada2_file  <- args[3]
output_dir  <- args[4]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Read taxonomy tables ───────────────────────────────────────────────────────
read_tax <- function(path, pipeline_label) {
  if (!file.exists(path)) {
    warning("File not found: ", path, " — treating as empty for ", pipeline_label)
    return(data.frame(lca_taxon = character(0), lca_rank = character(0),
                      pipeline = character(0)))
  }
  read.delim(path, stringsAsFactors = FALSE) |>
    mutate(pipeline = pipeline_label)
}

otu   <- read_tax(otu_file,   "NGSpeciesID_OTU")
asv   <- read_tax(asv_file,   "NGSpeciesID_ASV")
dada2 <- read_tax(dada2_file, "DADA2_ASV")

# ── Species detected by each pipeline ─────────────────────────────────────────
species_by_pipeline <- bind_rows(
  otu   |> filter(lca_rank == "species") |> distinct(lca_taxon, pipeline),
  asv   |> filter(lca_rank == "species") |> distinct(lca_taxon, pipeline),
  dada2 |> filter(lca_rank == "species") |> distinct(lca_taxon, pipeline)
)

# ── Build species × pipeline matrix ───────────────────────────────────────────
species_matrix <- species_by_pipeline |>
  mutate(detected = TRUE) |>
  pivot_wider(names_from  = pipeline,
              values_from = detected,
              values_fill = FALSE) |>
  # Ensure all three columns exist even if a pipeline produced no results
  {function(df) {
    for (col in c("NGSpeciesID_OTU", "NGSpeciesID_ASV", "DADA2_ASV")) {
      if (!col %in% names(df)) df[[col]] <- FALSE
    }
    df
  }}() |>
  mutate(
    concordant = NGSpeciesID_OTU & NGSpeciesID_ASV & DADA2_ASV,
    asv_only   = NGSpeciesID_ASV & !NGSpeciesID_OTU,
    dada2_only = DADA2_ASV       & !NGSpeciesID_ASV,
    otu_only   = NGSpeciesID_OTU & !NGSpeciesID_ASV
  ) |>
  arrange(desc(concordant), desc(NGSpeciesID_OTU & NGSpeciesID_ASV),
          lca_taxon)

# ── Genus-level comparisons ────────────────────────────────────────────────────
genera_matrix <- bind_rows(
  otu   |> filter(lca_rank == "genus") |> distinct(lca_taxon, pipeline),
  asv   |> filter(lca_rank == "genus") |> distinct(lca_taxon, pipeline),
  dada2 |> filter(lca_rank == "genus") |> distinct(lca_taxon, pipeline)
) |>
  mutate(detected = TRUE) |>
  pivot_wider(names_from  = pipeline,
              values_from = detected,
              values_fill = FALSE) |>
  {function(df) {
    for (col in c("NGSpeciesID_OTU", "NGSpeciesID_ASV", "DADA2_ASV")) {
      if (!col %in% names(df)) df[[col]] <- FALSE
    }
    df
  }}()

# ── Write outputs ──────────────────────────────────────────────────────────────
concordance_csv <- file.path(output_dir, "species_concordance.csv")
write.csv(species_matrix, concordance_csv, row.names = FALSE)
message("Species concordance written to: ", concordance_csv)

genera_csv <- file.path(output_dir, "genera_concordance.csv")
write.csv(genera_matrix, genera_csv, row.names = FALSE)
message("Genus concordance written to: ", genera_csv)

# ── Summary report ─────────────────────────────────────────────────────────────
summary_lines <- c(
  "=== ASV vs OTU SPECIES DETECTION SUMMARY ===",
  sprintf("Total species (union):              %d", nrow(species_matrix)),
  sprintf("Concordant (all 3 pipelines):       %d", sum(species_matrix$concordant)),
  sprintf("ASV-only (NGSpeciesID ASV > OTU):   %d", sum(species_matrix$asv_only)),
  sprintf("DADA2-only:                         %d", sum(species_matrix$dada2_only)),
  sprintf("OTU-only (lost in ASV transition):  %d", sum(species_matrix$otu_only)),
  "",
  "--- Concordant species ---",
  if (any(species_matrix$concordant))
    paste0("  ", species_matrix$lca_taxon[species_matrix$concordant])
  else
    "  (none)",
  "",
  "--- ASV-only species (investigate — genuine or artifact?) ---",
  if (any(species_matrix$asv_only))
    paste0("  ", species_matrix$lca_taxon[species_matrix$asv_only])
  else
    "  (none)",
  "",
  "--- DADA2-only species (verify manually via BLAST) ---",
  if (any(species_matrix$dada2_only))
    paste0("  ", species_matrix$lca_taxon[species_matrix$dada2_only])
  else
    "  (none)",
  "",
  "--- OTU-only species (likely OTU consensus blurring artifacts) ---",
  if (any(species_matrix$otu_only))
    paste0("  ", species_matrix$lca_taxon[species_matrix$otu_only])
  else
    "  (none)",
  "",
  "--- ASV-to-OTU mapping (manual step) ---",
  "To count NGSpeciesID ASVs per OTU cluster, run:",
  "  cd-hit-est -i asv_representatives.fasta -o asv_clustered_97.fasta -c 0.97",
  "  Then parse asv_clustered_97.fasta.clstr"
)

summary_path <- file.path(output_dir, "concordance_summary.txt")
writeLines(summary_lines, summary_path)
cat(paste(summary_lines, collapse = "\n"), "\n")
message("\nSummary written to: ", summary_path)
