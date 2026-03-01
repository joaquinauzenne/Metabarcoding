#' Initialize IUCN API key
#'
#' @param key character; IUCN Red List API token (defaults to IUCN_REDLIST_KEY env var)
#' @return invisible NULL
#' @export
rl_init <- function(key = Sys.getenv("IUCN_REDLIST_KEY")) {
  if (!nzchar(key)) {
    stop("Please set your IUCN API key via rl_init(key) or the IUCN_REDLIST_KEY env var.")
  }
  options(iucn_redlist_key = key)
  invisible(NULL)
}


#' Fetch species range shapefiles from IUCN
#'
#' @param species character vector; scientific names
#' @param path directory to download & unzip (defaults to tempdir())
#' @return named list of shapefile paths
#' @export
fetch_species_ranges <- function(species, path = tempdir()) {
  key <- getOption("iucn_redlist_key", Sys.getenv("IUCN_REDLIST_KEY"))
  if (!nzchar(key)) stop("IUCN API key not set. Run rl_init().")

  fs::dir_create(path)  # ensure directory exists :contentReference[oaicite:6]{index=6}
  out <- list()

  for (sp in species) {
    # Resolve taxon ID via rredlist
    info <- rredlist::rl_species(sp, key = key)$result
    if (nrow(info) == 0) stop("Species not found: ", sp)
    id <- info$taxonid[1]

    # Download zipped shapefile from Spatial API
    url <- sprintf("https://apiv3.iucnredlist.org/spatial/species/id/%s?token=%s", id, key)
    resp <- httr::GET(url)
    httr::stop_for_status(resp)

    zipfile <- fs::path(path, paste0(gsub(" ", "_", sp), ".zip"))
    writeBin(httr::content(resp, "raw"), zipfile)
    unzip_dir <- fs::path(path, gsub(" ", "_", sp))
    unzip(zipfile, exdir = unzip_dir)

    shp <- fs::dir_ls(unzip_dir, glob = "*.shp")[1]
    out[[sp]] <- shp
  }
  out
}


#' Read a range shapefile into sf
#'
#' @param file path to .shp
#' @return sf POLYGON/MULTIPOLYGON object
#' @export
read_range_shapefile <- function(file) {
  sf::st_read(file, quiet = TRUE)
}


#' Combine multiple species range sf objects
#'
#' @param sf_list named list of sf objects
#' @return single sf with `species` column
#' @export
combine_ranges <- function(sf_list) {
  purrr::imap_dfr(sf_list, ~{
    sf::st_set_geometry(.x, "geometry") %>%
      dplyr::mutate(species = .y)
  })
}


#' Crop ranges to a study boundary
#'
#' @param ranges_sf sf of species ranges
#' @param area_sf sf polygon for study area
#' @return cropped sf
#' @export
crop_to_study_area <- function(ranges_sf, area_sf) {
  sf::st_intersection(ranges_sf, sf::st_make_valid(area_sf))
}


#' Plot species ranges with leaflet
#'
#' @param ranges_sf sf with `species` and geometry
#' @param area_sf study boundary sf
#' @return leaflet map
#' @export
plot_range_map <- function(ranges_sf, area_sf) {
  pal <- leaflet::colorFactor("viridis", domain = ranges_sf$species)
  leaflet::leaflet() %>%
    leaflet::addTiles() %>%
    leaflet::addPolygons(data = area_sf, color = "blue", weight = 2, fill = FALSE) %>%
    leaflet::addPolygons(data = ranges_sf,
                         color = ~pal(species),
                         weight = 1,
                         opacity = 0.8,
                         group = ~species) %>%
    leaflet::addLayersControl(
      overlayGroups = unique(ranges_sf$species),
      options = leaflet::layersControlOptions(collapsed = FALSE)
    )
}


#' Get top N species by count for a barcode
#'
#' @param summary_df data.frame with `barcode`,`species`,`count`
#' @param barcode character
#' @param n integer
#' @return character vector
#' @export
get_top_species <- function(summary_df, barcode, n = 5) {
  summary_df %>%
    dplyr::filter(barcode == barcode) %>%
    dplyr::arrange(dplyr::desc(count)) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(species)
}


#' Load a study-area shapefile
#'
#' @param path to shapefile
#' @return sf polygon
#' @export
get_study_area <- function(path) {
  sf::st_read(path, quiet = TRUE)
}


#' Cache downloaded shapefiles to a directory
#'
#' @param file_paths character vector of shapefile paths
#' @param cache_dir directory to store cache
#' @return invisible NULL
#' @export
cache_range_data <- function(file_paths, cache_dir = "iucn_cache") {
  fs::dir_create(cache_dir)
  fs::file_copy(file_paths, cache_dir, overwrite = FALSE)
  invisible(NULL)
}


#' Clear cached shapefiles
#'
#' @export
clean_cache <- function(cache_dir = "iucn_cache") {
  if (fs::dir_exists(cache_dir)) fs::dir_delete(cache_dir)
}


