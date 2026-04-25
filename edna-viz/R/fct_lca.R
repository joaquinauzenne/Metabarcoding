#' In-app LCA (Lowest Common Ancestor) taxonomy assignment
#'
#' Mirrors the logic of `scripts/lca_assignment.R` as a pure R function
#' suitable for use inside the Shiny app's reactive pipeline.
#'
#' @param blast_df Canonical BLAST data.frame from [parse_blast_web()] or
#'   [parse_blast_cli()].
#' @param min_pident Minimum percent identity threshold (default 97).
#' @param min_qcov Minimum query coverage threshold (default 90).
#' @param bitscore_margin Keep hits within this fraction of the best bitscore
#'   (default 0.98, i.e., within 2%).
#' @return data.frame with one row per query, columns:
#'   query_id, lca_taxon, lca_rank, confidence,
#'   best_pident, best_qcov, n_species, all_species, common_names
#' @noRd
apply_lca <- function(blast_df,
                      min_pident      = 97,
                      min_qcov        = 90,
                      bitscore_margin = 0.98) {

  if (nrow(blast_df) == 0) {
    return(
      data.frame(query_id = character(0), lca_taxon = character(0),
                 lca_rank = character(0), confidence = character(0),
                 best_pident = numeric(0), best_qcov = numeric(0),
                 n_species = integer(0), all_species = character(0),
                 common_names = character(0), stringsAsFactors = FALSE)
    )
  }

  blast_filtered <- blast_df |>
    dplyr::filter(
      !is.na(pident), !is.na(qcov),
      pident >= min_pident,
      qcov   >= min_qcov
    ) |>
    # Split semicolon-delimited scientific names (can happen in CLI format)
    tidyr::separate_rows(scientific_name, sep = ";") |>
    dplyr::mutate(scientific_name = trimws(scientific_name)) |>
    dplyr::filter(
      !is.na(scientific_name),
      scientific_name != "",
      !grepl("uncultured|environmental|synthetic|unclassified",
             scientific_name, ignore.case = TRUE)
    )

  if (nrow(blast_filtered) == 0) {
    return(
      data.frame(query_id = character(0), lca_taxon = character(0),
                 lca_rank = character(0), confidence = character(0),
                 best_pident = numeric(0), best_qcov = numeric(0),
                 n_species = integer(0), all_species = character(0),
                 common_names = character(0), stringsAsFactors = FALSE)
    )
  }

  lca_results <- blast_filtered |>
    dplyr::group_by(query_id) |>
    dplyr::filter(
      is.na(bitscore) | bitscore >= max(bitscore, na.rm = TRUE) * bitscore_margin
    ) |>
    dplyr::summarize(
      n_top_hits    = dplyr::n(),
      n_species     = dplyr::n_distinct(scientific_name),
      best_pident   = max(pident, na.rm = TRUE),
      best_qcov     = max(qcov,   na.rm = TRUE),
      best_species  = if (dplyr::n_distinct(scientific_name) == 1L)
                        dplyr::first(scientific_name) else NA_character_,
      all_species   = paste(unique(scientific_name), collapse = "; "),
      genera        = paste(unique(sub(" .*", "", scientific_name)), collapse = "; "),
      n_genera      = dplyr::n_distinct(sub(" .*", "", scientific_name)),
      common_names  = paste(
        unique(common_name[!is.na(common_name) &
                             common_name != "" &
                             common_name != "N/A"]),
        collapse = "; "
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      lca_rank = dplyr::case_when(
        n_species == 1L ~ "species",
        n_genera  == 1L ~ "genus",
        TRUE            ~ "family_or_higher"
      ),
      lca_taxon = dplyr::case_when(
        lca_rank == "species" ~ best_species,
        lca_rank == "genus"   ~ sub(";.*", "", genera),
        TRUE                  ~ "unresolved"
      ),
      confidence = dplyr::case_when(
        best_pident >= 99 & lca_rank == "species"                 ~ "high",
        best_pident >= 97 & lca_rank %in% c("species", "genus")   ~ "medium",
        TRUE                                                       ~ "low"
      )
    ) |>
    dplyr::select(query_id, lca_taxon, lca_rank, confidence,
                  best_pident, best_qcov, n_top_hits, n_species, all_species, common_names)

  lca_results
}

#' Summarise LCA results for display
#'
#' Groups by lca_taxon and counts queries per species, for use in the filter
#' module bar chart and species selector.
#'
#' @param lca_df Output of [apply_lca()].
#' @param min_confidence Minimum confidence level to include. One of
#'   "high", "medium", "low". Default "low" (include all).
#' @return data.frame: scientific_name, n_queries, confidence_max
#' @noRd
summarise_lca <- function(lca_df, min_confidence = "low") {
  conf_levels <- c("low" = 1L, "medium" = 2L, "high" = 3L)
  min_level   <- conf_levels[[min_confidence]]

  lca_df |>
    dplyr::filter(
      !is.na(lca_taxon),
      lca_taxon != "unresolved",
      lca_rank  == "species",
      conf_levels[confidence] >= min_level
    ) |>
    dplyr::group_by(scientific_name = lca_taxon) |>
    dplyr::summarize(
      n_queries      = dplyr::n(),
      confidence_max = names(which.max(conf_levels[confidence])),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_queries))
}
