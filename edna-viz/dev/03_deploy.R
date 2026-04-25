## 03_deploy.R — Deployment commands
##
## Run these commands to build and deploy the edna-viz app.
## Do NOT source this file — run lines individually as needed.
## Working directory should be the edna-viz/ package root.

# ── Build the package ─────────────────────────────────────────────────────────

# Install the package locally (from the Metabarcoding/ repo root)
# devtools::install_local("edna-viz")

# Build a source tarball for distribution
# devtools::build()

# ── Docker ────────────────────────────────────────────────────────────────────

# The Dockerfile lives at edna-viz/Dockerfile (rocker/shiny:4.3.2 base).
# Build and run from the edna-viz/ directory:
#
#   docker build -t edna-viz .
#   docker run -p 3838:3838 edna-viz
#   # Open http://localhost:3838
#
# Or add the Golem Dockerfile scaffold (will overwrite existing Dockerfile):
# golem::add_dockerfile()

# ── shinyapps.io ─────────────────────────────────────────────────────────────

# Deploy to shinyapps.io (requires rsconnect configured with your account token)
# rsconnect::setAccountInfo(name = "<ACCOUNT>", token = "<TOKEN>", secret = "<SECRET>")
# rsconnect::deployApp(appDir = ".", appName = "edna-viz")

# ── Posit Connect ─────────────────────────────────────────────────────────────
# golem::add_dockerfile_with_renv()  # generates a renv-based Dockerfile
# rsconnect::deployApp()

# ── IUCN spatial data (one-time, before first deploy) ─────────────────────────
# Run the preprocessing script to generate mammal_ranges.rds:
#
#   Rscript data-raw/prep_iucn_ranges.R \
#     /path/to/TERRESTRIAL_MAMMALS.shp \
#     /path/to/taxonomy.tsv
#
# The .rds file is saved to inst/app/app_data/mammal_ranges.rds.
# It is excluded from git (*.rds in .gitignore) — regenerate on each machine.
