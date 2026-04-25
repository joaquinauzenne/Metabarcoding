library(testthat)

# Helper: minimal canonical BLAST data.frame for testing
make_blast_df <- function() {
  data.frame(
    query_id        = c("q1", "q1", "q2", "q3", "q4"),
    subject_id      = c("A1", "A2", "B1", "C1", "D1"),
    pident          = c(99.5, 99.2, 97.1, 95.0, 99.9),
    qcov            = c(99.0, 98.5, 91.0, 88.0, 99.0),
    evalue          = c(1e-60, 1e-59, 5e-45, 1e-30, 1e-70),
    bitscore        = c(250.0, 248.0, 200.0, 150.0, 300.0),
    taxid           = c(9503L,  9503L,  9600L,  1234L, 9999L),
    scientific_name = c("Alouatta seniculus", "Alouatta seniculus",
                        "Ateles hybridus",   "Dasypus novemcinctus",
                        "uncultured mammal"),
    common_name     = c("howler monkey", "howler monkey",
                        "brown spider monkey", "nine-banded armadillo",
                        "NA"),
    description     = c("desc1", "desc2", "desc3", "desc4", "desc5"),
    stringsAsFactors = FALSE
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Threshold filtering
# ─────────────────────────────────────────────────────────────────────────────
test_that("apply_lca retains rows above pident threshold", {
  df  <- make_blast_df()
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90)
  # q3 has pident=95 < 97 — should be excluded
  expect_false("Dasypus novemcinctus" %in% lca$lca_taxon)
})

test_that("apply_lca retains rows above qcov threshold", {
  df  <- make_blast_df()
  lca <- apply_lca(df, min_pident = 90, min_qcov = 90)
  # q3 has qcov=88 < 90 — should be excluded
  expect_false("Dasypus novemcinctus" %in% lca$lca_taxon)
})

test_that("apply_lca excludes uncultured/environmental hits", {
  df  <- make_blast_df()
  lca <- apply_lca(df, min_pident = 0, min_qcov = 0)
  expect_false("uncultured mammal" %in% lca$lca_taxon)
})

# ─────────────────────────────────────────────────────────────────────────────
# LCA rank assignment
# ─────────────────────────────────────────────────────────────────────────────
test_that("apply_lca: single-species top hits → species rank", {
  df  <- make_blast_df()
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90)

  # q1 has two hits both to Alouatta seniculus → species rank
  q1 <- lca[lca$query_id == "q1", ]
  expect_equal(nrow(q1), 1L)
  expect_equal(q1$lca_rank, "species")
  expect_equal(q1$lca_taxon, "Alouatta seniculus")
})

test_that("apply_lca: same genus but different species → genus rank", {
  # Two hits to different Alouatta species within bitscore margin
  df <- data.frame(
    query_id        = c("qX", "qX"),
    subject_id      = c("S1", "S2"),
    pident          = c(98.0, 97.5),
    qcov            = c(99.0, 98.0),
    evalue          = c(1e-60, 2e-59),
    bitscore        = c(250.0, 248.0),   # within 2% margin
    taxid           = c(9503L, 9600L),
    scientific_name = c("Alouatta seniculus", "Alouatta palliata"),
    common_name     = c("howler monkey", "mantled howler"),
    description     = c("d1", "d2"),
    stringsAsFactors = FALSE
  )
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90)
  expect_equal(lca$lca_rank, "genus")
  expect_equal(lca$lca_taxon, "Alouatta")
})

test_that("apply_lca: different genera → family_or_higher rank", {
  df <- data.frame(
    query_id        = c("qY", "qY"),
    subject_id      = c("S1", "S2"),
    pident          = c(98.0, 97.5),
    qcov            = c(99.0, 98.0),
    evalue          = c(1e-60, 2e-59),
    bitscore        = c(250.0, 249.0),
    taxid           = c(9503L, 1234L),
    scientific_name = c("Alouatta seniculus", "Dasypus novemcinctus"),
    common_name     = c("howler monkey", "armadillo"),
    description     = c("d1", "d2"),
    stringsAsFactors = FALSE
  )
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90)
  expect_equal(lca$lca_rank, "family_or_higher")
  expect_equal(lca$lca_taxon, "unresolved")
})

# ─────────────────────────────────────────────────────────────────────────────
# Confidence scoring
# ─────────────────────────────────────────────────────────────────────────────
test_that("apply_lca: pident ≥ 99 + species rank → high confidence", {
  df  <- make_blast_df()
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90)
  q1  <- lca[lca$query_id == "q1", ]
  expect_equal(q1$confidence, "high")
})

test_that("apply_lca: 97 ≤ pident < 99 + species rank → medium confidence", {
  df  <- make_blast_df()
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90)
  q2  <- lca[lca$query_id == "q2", ]
  expect_equal(q2$confidence, "medium")
})

# ─────────────────────────────────────────────────────────────────────────────
# summarise_lca
# ─────────────────────────────────────────────────────────────────────────────
test_that("summarise_lca returns one row per species", {
  df  <- make_blast_df()
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90)
  ss  <- summarise_lca(lca)
  expect_true(all(!duplicated(ss$scientific_name)))
})

test_that("summarise_lca: confidence filter works", {
  df  <- make_blast_df()
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90)
  ss_all  <- summarise_lca(lca, min_confidence = "low")
  ss_high <- summarise_lca(lca, min_confidence = "high")
  expect_true(nrow(ss_high) <= nrow(ss_all))
})

# ─────────────────────────────────────────────────────────────────────────────
# Edge cases
# ─────────────────────────────────────────────────────────────────────────────
test_that("apply_lca handles empty input gracefully", {
  df  <- make_blast_df()[0, ]
  lca <- apply_lca(df)
  expect_equal(nrow(lca), 0L)
  expect_true(all(CANONICAL_COLS[1:2] %in% c("query_id", "subject_id")))
})

test_that("apply_lca: bitscore_margin keeps only top hits", {
  # q1 has bitscore 250 and 248 — both within 2% margin (250*0.98 = 245)
  df <- make_blast_df()
  lca <- apply_lca(df, min_pident = 97, min_qcov = 90, bitscore_margin = 0.98)
  q1 <- lca[lca$query_id == "q1", ]
  expect_equal(q1$n_top_hits, 2L)
})
