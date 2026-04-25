## 02_dev.R — Development cycle commands
##
## Run these commands during active development of edna-viz.
## Do NOT source this file — run lines individually as needed.
## Working directory should be the edna-viz/ package root.

# ── Load / reload ─────────────────────────────────────────────────────────────

# Load all package functions without installing (fast, use during development)
devtools::load_all()

# ── Documentation ─────────────────────────────────────────────────────────────

# Regenerate NAMESPACE and man/ from roxygen2 comments
devtools::document()

# ── Testing ───────────────────────────────────────────────────────────────────

# Run all tests
devtools::test()

# Run tests with detailed reporter
devtools::test(reporter = testthat::ProgressReporter)

# Run a specific test file
# testthat::test_file("tests/testthat/test-parse_blast.R")
# testthat::test_file("tests/testthat/test-filter_logic.R")

# ── Checking ──────────────────────────────────────────────────────────────────

# Full R CMD check (use before release or PR)
devtools::check()

# Quick check without building vignettes/manual
devtools::check(document = FALSE, vignettes = FALSE)

# ── Run the app ───────────────────────────────────────────────────────────────

# Launch the app in development mode (hot-reloads on file save if using golem)
golem::run_dev()

# Or launch directly via shiny
# shiny::runApp(".")

# ── Spell check ───────────────────────────────────────────────────────────────
# spelling::spell_check_package()

# ── Code coverage ─────────────────────────────────────────────────────────────
# covr::package_coverage()
# covr::report()
