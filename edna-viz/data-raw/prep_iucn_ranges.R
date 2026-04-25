# prep_iucn_ranges.R — One-time preprocessing of IUCN Terrestrial Mammals shapefile
#
# Run this script ONCE locally after downloading the IUCN spatial data.
# Output: edna-viz/inst/app/app_data/mammal_ranges.rds
#
# Usage:
#   # With hardcoded seed list only:
#   Rscript data-raw/prep_iucn_ranges.R /path/to/TERRESTRIAL_MAMMALS.shp
#
#   # With species derived from BLAST results (recommended — run after each pipeline job):
#   Rscript data-raw/prep_iucn_ranges.R /path/to/TERRESTRIAL_MAMMALS.shp \
#           eDNA/consensus/combined/taxonomy/taxonomy.tsv
#
# When taxonomy.tsv is provided, all species-rank LCA assignments are added to
# the target list automatically. The seed list below is always included as a
# fallback so focal species appear in the map even before pipeline results exist.
#
# Prerequisites:
#   1. Register at https://www.iucnredlist.org/resources/spatial-data-download
#   2. Download "Terrestrial Mammals" spatial data (~594 MB zipped)
#   3. Unzip to a local directory (not tracked by git)
#
# The script filters to target species, simplifies geometries for performance,
# and saves a compact .rds for bundling with the app.

library(sf)
library(rmapshaper)
library(dplyr)

# ── Command-line arguments ─────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  # Interactive fallback: use a hardcoded default path
  SHAPEFILE_PATH <- "~/Downloads/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp"
  TAXONOMY_TSV   <- NULL
  message("No arguments provided — using default SHAPEFILE_PATH.")
  message("Usage: Rscript prep_iucn_ranges.R <shapefile.shp> [taxonomy.tsv]")
} else {
  SHAPEFILE_PATH <- args[1]
  TAXONOMY_TSV   <- if (length(args) >= 2) args[2] else NULL
}

# ── Seed species list ──────────────────────────────────────────────────────────
# Always included regardless of taxonomy.tsv. Focal species for the study area
# (Tatacoa Desert and Magdalena River Valley, Colombia).
SEED_MAMMALS <- c(
  "Alouatta seniculus",        # Colombian red howler monkey — confirmed in field
  "Ateles hybridus",           # Brown spider monkey (Critically Endangered)
  "Cebus versicolor",          # Varied white-fronted capuchin
  "Dasypus novemcinctus",      # Nine-banded armadillo
  "Didelphis marsupialis",     # Common opossum
  "Tamandua mexicana",         # Northern tamandua
  "Mazama americana",          # Red brocket deer
  "Pecari tajacu",             # Collared peccary
  "Nasua nasua",               # South American coati
  "Cerdocyon thous",           # Crab-eating fox
  "Leopardus pardalis",        # Ocelot
  "Tapirus terrestris"         # South American tapir
)

# ── Extend target list from BLAST results ─────────────────────────────────────
blast_species <- character(0)
if (!is.null(TAXONOMY_TSV)) {
  if (!file.exists(TAXONOMY_TSV)) {
    warning("taxonomy.tsv not found at: ", TAXONOMY_TSV, " — using seed list only.")
  } else {
    tax <- read.delim(TAXONOMY_TSV, stringsAsFactors = FALSE)
    if (!"lca_taxon" %in% names(tax) || !"lca_rank" %in% names(tax)) {
      warning("taxonomy.tsv missing 'lca_taxon' or 'lca_rank' columns — using seed list only.")
    } else {
      blast_species <- tax |>
        dplyr::filter(lca_rank == "species", !is.na(lca_taxon)) |>
        dplyr::pull(lca_taxon) |>
        unique()
      message(sprintf(
        "Loaded %d species-rank detections from taxonomy.tsv: %s",
        length(blast_species),
        paste(blast_species, collapse = ", ")
      ))
    }
  }
}

TARGET_MAMMALS <- unique(c(SEED_MAMMALS, blast_species))
message(sprintf("Total target species: %d", length(TARGET_MAMMALS)))

# ── Load and filter IUCN shapefile ────────────────────────────────────────────
if (!file.exists(SHAPEFILE_PATH)) {
  stop("Shapefile not found at: ", SHAPEFILE_PATH,
       "\nUpdate SHAPEFILE_PATH in this script.")
}

message("Reading IUCN terrestrial mammals shapefile...")
mammals_raw <- sf::st_read(SHAPEFILE_PATH, quiet = TRUE)

# The sci_name field was renamed from "binomial" in 2023 IUCN data release.
# Detect which column is used.
name_col <- if ("sci_name" %in% names(mammals_raw)) "sci_name" else "binomial"
message(sprintf("Using species name column: '%s'", name_col))

message(sprintf("Total polygons in shapefile: %d", nrow(mammals_raw)))

# Filter:
#   presence == 1  → Extant
#   origin   == 1  → Native
#   seasonal %in% c(1, 2) → Resident or Breeding season
mammals_filtered <- mammals_raw |>
  dplyr::filter(
    .data[[name_col]] %in% TARGET_MAMMALS,
    presence %in% 1L,
    origin   %in% 1L,
    seasonal %in% c(1L, 2L)
  ) |>
  sf::st_transform(crs = 4326)

n_found    <- length(unique(mammals_filtered[[name_col]]))
n_target   <- length(TARGET_MAMMALS)
n_missing  <- TARGET_MAMMALS[!TARGET_MAMMALS %in% mammals_filtered[[name_col]]]

message(sprintf("Found %d / %d target species in shapefile.", n_found, n_target))
if (length(n_missing) > 0) {
  message("Missing from shapefile (no IUCN range data):")
  message(paste0("  - ", n_missing, collapse = "\n"))
}

# ── Simplify geometries for web performance ───────────────────────────────────
# keep = 0.05 retains 5% of vertices — sufficient for country/region scale.
# keep_shapes = TRUE prevents polygons from disappearing at low vertex counts.
message("Simplifying geometries (keep = 0.05)...")
mammals_simplified <- rmapshaper::ms_simplify(
  mammals_filtered,
  keep        = 0.05,
  keep_shapes = TRUE
)

# ── Save output ────────────────────────────────────────────────────────────────
out_path <- file.path(
  dirname(dirname(rstudioapi::getActiveDocumentContext()$path)),
  "inst", "app", "app_data", "mammal_ranges.rds"
)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(mammals_simplified, out_path)

size_kb <- file.size(out_path) / 1024
message(sprintf(
  "Saved %d polygons for %d species (%.1f KB) to:\n  %s",
  nrow(mammals_simplified), n_found, size_kb, out_path
))
message("\nDone. Re-run this script after adding species to TARGET_MAMMALS.")
