#' Fetch GBIF occurrence records for a set of species
#'
#' Queries the GBIF Occurrence API via [rgbif::occ_search()] for each species
#' name in `species_names`. Returns a combined data.frame suitable for
#' rendering as Leaflet circle markers.
#'
#' @param species_names Character vector of scientific names.
#' @param country ISO2 country code to restrict results (default "CO" for
#'   Colombia — Tatacoa / Magdalena Valley study area).
#' @param limit Maximum records per species (default 200).
#' @return data.frame with columns: scientific_name, decimalLatitude,
#'   decimalLongitude, year, datasetName, occurrenceID. Rows with missing
#'   coordinates are dropped.
#' @noRd
fetch_gbif_occurrences <- function(species_names,
                                   country = "CO",
                                   limit   = 200L) {
  if (length(species_names) == 0) return(empty_occurrences())

  results <- lapply(species_names, function(sp) {
    tryCatch({
      res <- rgbif::occ_search(
        scientificName = sp,
        country        = country,
        hasCoordinate  = TRUE,
        limit          = limit,
        fields         = c("scientificName", "decimalLatitude", "decimalLongitude",
                           "year", "datasetName", "occurrenceID")
      )
      if (is.null(res$data) || nrow(res$data) == 0L) return(NULL)
      df <- res$data
      df$species_query <- sp
      df
    }, error = function(e) {
      warning("GBIF query failed for '", sp, "': ", conditionMessage(e))
      NULL
    })
  })

  combined <- dplyr::bind_rows(Filter(Negate(is.null), results))
  if (nrow(combined) == 0L) return(empty_occurrences())

  combined |>
    dplyr::filter(
      !is.na(decimalLatitude),
      !is.na(decimalLongitude)
    ) |>
    dplyr::rename(scientific_name = species_query) |>
    dplyr::select(scientific_name, decimalLatitude, decimalLongitude,
                  dplyr::any_of(c("year", "datasetName", "occurrenceID")))
}

#' Build a GBIF Maps API tile URL for a species
#'
#' Returns a lightweight tile URL for use with [leaflet::addTiles()] instead
#' of downloading individual occurrence points. Useful when species has
#' hundreds of thousands of records.
#'
#' @param taxon_key GBIF taxon key (integer). Obtain via
#'   [rgbif::name_backbone(name = "Species name")$usageKey].
#' @param style Tile style (default "purpleYellow.point").
#' @return URL string suitable for [leaflet::addTiles(urlTemplate = ...)].
#' @noRd
gbif_tile_url <- function(taxon_key,
                          style = "purpleYellow.point") {
  paste0(
    "https://api.gbif.org/v2/map/occurrence/density/{z}/{x}/{y}@1x.png?",
    "style=", style,
    "&taxonKey=", taxon_key
  )
}

#' @noRd
empty_occurrences <- function() {
  data.frame(
    scientific_name    = character(0),
    decimalLatitude    = numeric(0),
    decimalLongitude   = numeric(0),
    stringsAsFactors   = FALSE
  )
}
