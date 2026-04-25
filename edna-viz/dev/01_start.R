## 01_start.R — One-time package setup commands
##
## Run these commands manually once when setting up the edna-viz package
## for the first time or adding new modules/functions.
## Do NOT source this file — run lines individually as needed.

# Set the working directory to the edna-viz package root
# setwd("edna-viz")  # or open edna-viz.Rproj

# ── Package infrastructure ────────────────────────────────────────────────────

# Initialise Golem options (run once if not already in app_config.R)
# golem::fill_desc(
#   pkg_name       = "ednaviz",
#   pkg_title      = "eDNA Metabarcoding Visualization App",
#   pkg_description = "Interactive Leaflet map and BLAST result browser for eDNA metabarcoding data.",
#   author_first_name  = "Joaquin",
#   author_last_name   = "Auzenne",
#   author_email       = "",
#   repo_url           = "https://github.com/joaquinauzenne/Metabarcoding"
# )

# Enable magrittr pipe (adds utils-pipe.R)
# usethis::use_pipe()

# ── Add modules ──────────────────────────────────────────────────────────────
# Use these commands to scaffold new Shiny modules:
# golem::add_module(name = "upload")   # creates R/mod_upload.R
# golem::add_module(name = "filter")   # creates R/mod_filter.R
# golem::add_module(name = "map")      # creates R/mod_map.R

# ── Add utility functions ─────────────────────────────────────────────────────
# Use these commands to scaffold new fct_*.R files:
# golem::add_fct("detect_format")     # creates R/fct_detect_format.R
# golem::add_fct("parse_blast")       # creates R/fct_parse_blast.R
# golem::add_fct("lca")               # creates R/fct_lca.R
# golem::add_fct("gbif_query")        # creates R/fct_gbif_query.R

# ── Package dependencies ──────────────────────────────────────────────────────
# Add new dependencies with:
# usethis::use_package("packagename")        # Imports
# usethis::use_package("packagename", type = "Suggests")

# ── Data preparation ──────────────────────────────────────────────────────────
# Add data-raw scripts:
# usethis::use_data_raw("prep_iucn_ranges")  # creates data-raw/prep_iucn_ranges.R

# ── Testing ───────────────────────────────────────────────────────────────────
# usethis::use_testthat()             # initialise testthat (already done)
# usethis::use_test("parse_blast")    # creates tests/testthat/test-parse_blast.R
# usethis::use_test("filter_logic")   # creates tests/testthat/test-filter_logic.R
