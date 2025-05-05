#' Initialize IUCN API key
#'
#' @param key character; your IUCN Red List API token
#' @return invisible NULL (sets option or env var)
#' @export
rl_init <- function(key = Sys.getenv("IUCN_REDLIST_KEY")) {
  options(iucn_redlist_key = key)
}

#' Fetch species range shapefiles
#'
#' @param species character vector; species scientific names or IDs
#' @param path     character; directory to save downloaded files
#' @return named list of file paths
#' @export
fetch_species_ranges <- function(species, path = tempdir()) { ... }

#' Read a downloaded range shapefile
#'
#' @param file path to a shapefile (.shp)
#' @return sf object of type POLYGON/MULTIPOLYGON
#' @export
read_range_shapefile <- function(file) { ... }

#' Combine multiple species sf objects
#'
#' @param sf_list named list of sf objects
#' @return single sf data frame with a `species` column
#' @export
combine_ranges <- function(sf_list) { ... }

#' Crop ranges to a study area
#'
#' @param ranges_sf sf data frame of species ranges
#' @param area      sf object defining study boundary
#' @return cropped sf object
#' @export
crop_to_study_area <- function(ranges_sf, area) { ... }

#' Plot species ranges over a basemap
#'
#' @param ranges_sf sf data frame with species ranges
#' @param area      sf object for boundary overlay
#' @return leaflet map widget
#' @export
plot_range_map <- function(ranges_sf, area) { ... }

#' Get top N species from a summary table
#'
#' @param summary_df data frame with columns `barcode`,`species`,`count`
#' @param barcode    character; barcode identifier
#' @param n          integer; number of top species
#' @return character vector of species
#' @export
get_top_species <- function(summary_df, barcode, n = 5) { ... }

#' Load a study-area shapefile
#'
#' @param path file path to shapefile
#' @return sf polygon
#' @export
get_study_area <- function(path) { ... }

#' Cache downloaded shapefiles
#'
#' @param file_paths character vector
#' @param cache_dir  directory to store
#' @return invisible cache index
#' @export
cache_range_data <- function(file_paths, cache_dir = "iucn_cache") { ... }

#' Clear cached shapefiles
#'
#' @export
clean_cache <- function() { ... }

