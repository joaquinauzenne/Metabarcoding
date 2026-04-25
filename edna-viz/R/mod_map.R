#' Map Module — Leaflet Map with IUCN Ranges, GBIF Points & eDNA Sites
#'
#' Renders an interactive Leaflet map with three overlay layers:
#' \describe{
#'   \item{Layer 1 (IUCN)}{Pre-processed IUCN range polygons from
#'     \code{inst/app/app_data/mammal_ranges.rds}}
#'   \item{Layer 2 (GBIF)}{Occurrence points fetched live from GBIF API,
#'     clustered for performance}
#'   \item{Layer 3 (eDNA sites)}{Hardcoded sampling site coordinates for the
#'     Tatacoa Desert and Magdalena River Valley study area}
#' }
#'
#' @param id Shiny module ID.
#' @param species_data Reactive returning the list from [mod_filter_server()].

# ── Sampling site coordinates (Magdalena River Valley, Colombia) ──
# Six field sites from the manuscript (DMS converted to decimal degrees):
#   1. Tatacoa Desert (Huila) — 3°12'00"N 75°02'48"W
#   2. Hacienda la Germania (Norcasia, Caldas) — 5°37'37"N 74°49'27"W
#   3. Hacienda Puerto Triunfo (Honda, Tolima) — 5°09'23"N 74°48'16"W
#   4. Hacienda San Juan de Carare (Cimitarra, Santander) — 6°41'44"N 74°09'24"W
#   5. Hacienda Lucitania (Cimitarra, Santander) — 6°26'15"N 74°07'12"W
#   6. Hacienda Antares (Puerto Parra, Santander) — 6°34'40"N 73°51'08"W
EDNA_SITES <- data.frame(
  site_name  = c(
    "Tatacoa Desert",
    "Hacienda la Germania",
    "Hacienda Puerto Triunfo",
    "Hacienda San Juan de Carare",
    "Hacienda Lucitania",
    "Hacienda Antares"
  ),
  latitude   = c(3.2000, 5.6269, 5.1564, 6.6956, 6.4375, 6.5778),
  longitude  = c(-75.0467, -74.8242, -74.8044, -74.1567, -74.1200, -73.8522),
  stringsAsFactors = FALSE
)

# ── UI ────────────────────────────────────────────────────────────────────────
mod_map_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    leaflet::leafletOutput(ns("map"), height = "600px"),
    shiny::div(
      class = "mt-2",
      shiny::checkboxInput(ns("show_gbif"),  "Show GBIF occurrences", value = TRUE),
      shiny::checkboxInput(ns("show_sites"), "Show eDNA sampling sites", value = TRUE)
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────
mod_map_server <- function(id, species_data) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Load pre-processed IUCN ranges once at startup
    mammal_ranges <- shiny::reactive({
      rds_path <- system.file(
        "app/app_data/mammal_ranges.rds",
        package = "ednaviz"
      )
      if (!file.exists(rds_path)) return(NULL)
      readRDS(rds_path)
    })

    # Species colour palette (shared across layers)
    species_pal <- shiny::reactive({
      sel <- species_data()$selected_species
      if (length(sel) == 0L) return(NULL)
      leaflet::colorFactor(
        palette = "Dark2",
        domain  = sel
      )
    })

    # Base map
    output$map <- leaflet::renderLeaflet({
      leaflet::leaflet() |>
        leaflet::addProviderTiles("CartoDB.Positron",
                                  group = "CartoDB Positron") |>
        leaflet::addProviderTiles("Esri.WorldImagery",
                                  group = "Esri World Imagery") |>
        leaflet::addProviderTiles("Esri.WorldTopoMap",
                                  group = "Esri Topo") |>
        leaflet::setView(lng = -74.5, lat = 5.0, zoom = 7) |>
        leaflet::addLayersControl(
          baseGroups    = c("CartoDB Positron", "Esri World Imagery", "Esri Topo"),
          overlayGroups = c("IUCN Ranges", "GBIF Occurrences", "eDNA Sites"),
          options       = leaflet::layersControlOptions(collapsed = FALSE)
        )
    })

    # ── Layer 1: IUCN range polygons ─────────────────────────────────────────
    shiny::observe({
      sel    <- species_data()$selected_species
      ranges <- mammal_ranges()
      pal    <- species_pal()

      leaflet::leafletProxy("map", session) |>
        leaflet::clearGroup("IUCN Ranges")

      if (is.null(ranges) || length(sel) == 0L) return()

      # Filter to selected species (IUCN shapefile uses sci_name column)
      name_col <- if ("sci_name" %in% names(ranges)) "sci_name" else "binomial"
      visible  <- ranges[ranges[[name_col]] %in% sel, ]

      if (nrow(visible) == 0L) return()

      leaflet::leafletProxy("map", session) |>
        leaflet::addPolygons(
          data        = visible,
          group       = "IUCN Ranges",
          fillColor   = ~pal(visible[[name_col]]),
          fillOpacity = 0.3,
          color       = ~pal(visible[[name_col]]),
          weight      = 1,
          popup       = ~paste0("<b>", visible[[name_col]], "</b><br>IUCN range")
        )
    })

    # ── Layer 2: GBIF occurrences (fetched live) ──────────────────────────────
    gbif_data <- shiny::eventReactive(
      list(species_data()$selected_species, input$show_gbif), {
        shiny::req(input$show_gbif)
        sel <- species_data()$selected_species
        if (length(sel) == 0L) return(NULL)
        fetch_gbif_occurrences(sel, country = "CO", limit = 200L)
      },
      ignoreNULL = FALSE
    )

    shiny::observe({
      occ <- gbif_data()
      pal <- species_pal()

      leaflet::leafletProxy("map", session) |>
        leaflet::clearGroup("GBIF Occurrences")

      if (!input$show_gbif || is.null(occ) || nrow(occ) == 0L || is.null(pal))
        return()

      leaflet::leafletProxy("map", session) |>
        leaflet::addCircleMarkers(
          data         = occ,
          lng          = ~decimalLongitude,
          lat          = ~decimalLatitude,
          group        = "GBIF Occurrences",
          color        = ~pal(scientific_name),
          radius       = 5,
          fillOpacity  = 0.7,
          stroke       = FALSE,
          popup        = ~paste0(
            "<b>", scientific_name, "</b><br>",
            "Year: ", ifelse(is.null(occ$year), "?", occ$year), "<br>",
            "Source: GBIF"
          ),
          clusterOptions = leaflet.extras::markerClusterOptions()
        )
    })

    # ── Layer 3: eDNA sampling sites ──────────────────────────────────────────
    shiny::observe({
      leaflet::leafletProxy("map", session) |>
        leaflet::clearGroup("eDNA Sites")

      if (!input$show_sites) return()

      leaflet::leafletProxy("map", session) |>
        leaflet::addCircleMarkers(
          data        = EDNA_SITES,
          lng         = ~longitude,
          lat         = ~latitude,
          group       = "eDNA Sites",
          color       = "#E74C3C",
          radius      = 8,
          fillOpacity = 0.9,
          stroke      = TRUE,
          weight      = 2,
          popup       = ~paste0("<b>", site_name, "</b><br>eDNA sampling site")
        )
    })
  })
}
