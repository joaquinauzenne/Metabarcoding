library(testthat)

# Helper: path to the bundled example file (web format)
demo_csv_path <- function() {
  # Try package inst/ first, then repo examples/
  p <- system.file("app/app_data/demo_blast.csv", package = "ednaviz")
  if (!file.exists(p)) {
    p <- file.path(
      system.file(package = "ednaviz"),
      "..", "..", "examples",
      "2ER9WHC1013-Alignment-Descriptions.csv"
    )
  }
  # Last resort: relative to test file
  if (!file.exists(p)) {
    p <- file.path(
      dirname(dirname(dirname(testthat::test_path()))),
      "examples",
      "2ER9WHC1013-Alignment-Descriptions.csv"
    )
  }
  normalizePath(p, mustWork = FALSE)
}

# Helper: write a temporary CLI TSV fixture
write_cli_fixture <- function(rows = NULL) {
  if (is.null(rows)) {
    rows <- c(
      paste(c("q1", "OM329056.1", "99.24", "99", "1e-59", "239",
              "9503", "Alouatta seniculus", "howler monkey",
              "Alouatta seniculus mitochondrion"), collapse = "\t"),
      paste(c("q2", "OM329046.1", "97.50", "92", "5e-55", "210",
              "9503;9600", "Alouatta seniculus;Alouatta palliata", "howler monkey",
              "Alouatta sp. mitochondrion"), collapse = "\t"),
      paste(c("q3", "AB123456.1", "85.00", "75", "1e-30", "150",
              "9000", "uncultured mammal", "NA",
              "uncultured environmental sample"), collapse = "\t")
    )
  }
  tmp <- tempfile(fileext = ".tsv")
  writeLines(rows, tmp)
  tmp
}

# ─────────────────────────────────────────────────────────────────────────────
# Format detection
# ─────────────────────────────────────────────────────────────────────────────
test_that("web format is detected from header containing 'Per. ident'", {
  csv <- demo_csv_path()
  skip_if_not(file.exists(csv), "Demo CSV not found")

  result <- detect_and_parse_blast(csv)
  expect_equal(result$format, "web")
})

test_that("CLI format is detected from headerless TSV", {
  tsv <- write_cli_fixture()
  on.exit(unlink(tsv))

  result <- detect_and_parse_blast(tsv)
  expect_equal(result$format, "cli")
})

# ─────────────────────────────────────────────────────────────────────────────
# Web format parsing
# ─────────────────────────────────────────────────────────────────────────────
test_that("web format: canonical columns are present", {
  csv <- demo_csv_path()
  skip_if_not(file.exists(csv), "Demo CSV not found")

  df <- parse_blast_web(csv)
  expect_true(all(CANONICAL_COLS %in% names(df)))
})

test_that("web format: percent signs are stripped from pident and qcov", {
  csv <- demo_csv_path()
  skip_if_not(file.exists(csv), "Demo CSV not found")

  df <- parse_blast_web(csv)
  expect_true(is.numeric(df$pident))
  expect_true(is.numeric(df$qcov))
  # 99% should become 99, not NA
  expect_true(all(df$pident >= 0 & df$pident <= 100, na.rm = TRUE))
  expect_true(all(df$qcov   >= 0 & df$qcov   <= 100, na.rm = TRUE))
})

test_that("web format: Alouatta seniculus is in scientific_name", {
  csv <- demo_csv_path()
  skip_if_not(file.exists(csv), "Demo CSV not found")

  df <- parse_blast_web(csv)
  expect_true("Alouatta seniculus" %in% df$scientific_name)
})

test_that("web format: no rows with empty scientific_name after parsing", {
  csv <- demo_csv_path()
  skip_if_not(file.exists(csv), "Demo CSV not found")

  df <- parse_blast_web(csv)
  expect_true(all(df$scientific_name != "" & !is.na(df$scientific_name)))
})

# ─────────────────────────────────────────────────────────────────────────────
# CLI format parsing
# ─────────────────────────────────────────────────────────────────────────────
test_that("CLI format: canonical columns are present", {
  tsv <- write_cli_fixture()
  on.exit(unlink(tsv))

  df <- parse_blast_cli(tsv)
  expect_true(all(CANONICAL_COLS %in% names(df)))
})

test_that("CLI format: semicolon-delimited taxids take the first entry", {
  tsv <- write_cli_fixture()
  on.exit(unlink(tsv))

  df <- parse_blast_cli(tsv)
  # Row q2 has "9503;9600" — should yield 9503
  q2 <- df[df$query_id == "q2", ]
  expect_equal(q2$taxid, 9503L)
  expect_equal(q2$scientific_name, "Alouatta seniculus")
})

test_that("CLI format: pident and qcov are numeric", {
  tsv <- write_cli_fixture()
  on.exit(unlink(tsv))

  df <- parse_blast_cli(tsv)
  expect_true(is.numeric(df$pident))
  expect_true(is.numeric(df$qcov))
})

test_that("CLI format: scientific notation e-values are parsed correctly", {
  rows <- paste(c("q1", "X1", "99.0", "99", "1e-120", "500",
                  "9503", "Alouatta seniculus", "howler monkey", "desc"),
                collapse = "\t")
  tsv <- write_cli_fixture(rows)
  on.exit(unlink(tsv))

  df <- parse_blast_cli(tsv)
  expect_equal(df$evalue, 1e-120)
})

# ─────────────────────────────────────────────────────────────────────────────
# Edge cases
# ─────────────────────────────────────────────────────────────────────────────
test_that("CLI format: uncultured/environmental rows are not dropped by parser", {
  tsv <- write_cli_fixture()
  on.exit(unlink(tsv))

  # The parser itself does NOT drop uncultured rows — LCA does
  df <- parse_blast_cli(tsv)
  expect_true("uncultured mammal" %in% df$scientific_name)
})

test_that("validate_blast_data returns empty vector for valid data", {
  csv <- demo_csv_path()
  skip_if_not(file.exists(csv), "Demo CSV not found")

  df <- parse_blast_web(csv)
  msgs <- validate_blast_data(df)
  expect_equal(length(msgs), 0L)
})

test_that("validate_blast_data flags missing canonical columns", {
  df <- data.frame(query_id = "q1", stringsAsFactors = FALSE)
  msgs <- validate_blast_data(df)
  expect_true(length(msgs) > 0L)
  expect_true(any(grepl("Missing columns", msgs)))
})
