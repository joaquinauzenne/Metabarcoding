#' Filter Module — Threshold Sliders, LCA Assignment & Bar Chart
#'
#' Consumes the canonical BLAST data.frame from [mod_upload_server()].
#' Applies LCA assignment, then lets the user filter by query coverage,
#' percent identity, and species selection.
#'
#' Hits below thresholds are flagged (reduced opacity in bar chart) rather
#' than dropped — important for 12S/16S reference database gaps in
#' Neotropical mammals.
#'
#' @param id Shiny module ID.
#' @param blast_data Reactive returning canonical BLAST data.frame (from
#'   [mod_upload_server()]).
#' @return (server) Reactive returning a list:
#'   \describe{
#'     \item{lca}{LCA-assigned data.frame}
#'     \item{selected_species}{Character vector of selected species names}
#'     \item{thresholds}{list(min_qcov, min_pident)}
#'   }

# ── UI ────────────────────────────────────────────────────────────────────────
mod_filter_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::hr(),
    shiny::h5("Filter thresholds"),
    shiny::sliderInput(
      ns("min_qcov"),
      label = "Min. Query Coverage (%)",
      min = 50, max = 100, value = 80, step = 1
    ),
    shiny::sliderInput(
      ns("min_pident"),
      label = "Min. % Identity",
      min = 85, max = 100, value = 95, step = 0.5
    ),
    shiny::hr(),
    shiny::h5("Species selection"),
    shiny::selectizeInput(
      ns("species_select"),
      label   = NULL,
      choices = NULL,
      multiple = TRUE,
      options = list(placeholder = "All species (or choose subset)")
    ),
    shiny::hr(),
    plotly::plotlyOutput(ns("bar_chart"), height = "300px")
  )
}

# ── Server ────────────────────────────────────────────────────────────────────
mod_filter_server <- function(id, blast_data) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Apply LCA to the uploaded data (all hits, before threshold filtering)
    lca_all <- shiny::reactive({
      bd <- blast_data()
      shiny::req(!is.null(bd), nrow(bd) > 0L)
      apply_lca(bd, min_pident = 0, min_qcov = 0)
    })

    # Summary per species (for bar chart and species selector)
    species_summary <- shiny::reactive({
      la <- lca_all()
      shiny::req(!is.null(la))
      # Flag rows that pass the current thresholds
      la |>
        dplyr::mutate(
          passes = best_pident >= input$min_pident & best_qcov >= input$min_qcov
        ) |>
        dplyr::filter(!is.na(lca_taxon), lca_taxon != "unresolved",
                      lca_rank == "species") |>
        dplyr::group_by(scientific_name = lca_taxon) |>
        dplyr::summarize(
          n_queries     = dplyr::n(),
          n_passing     = sum(passes),
          confidence    = names(which.max(c(low=1L, medium=2L, high=3L)[confidence])),
          .groups       = "drop"
        ) |>
        dplyr::arrange(dplyr::desc(n_passing), dplyr::desc(n_queries))
    })

    # Populate species selector
    shiny::observe({
      ss <- species_summary()
      shiny::req(!is.null(ss))
      choices <- stats::setNames(ss$scientific_name, ss$scientific_name)
      shiny::updateSelectizeInput(
        session, "species_select",
        choices  = choices,
        selected = character(0),
        server   = TRUE
      )
    })

    # Bar chart: detection frequency with above/below threshold colouring
    output$bar_chart <- plotly::renderPlotly({
      ss <- species_summary()
      shiny::req(!is.null(ss), nrow(ss) > 0L)

      # If user has selected species, filter to those
      if (length(input$species_select) > 0L) {
        ss <- ss |> dplyr::filter(scientific_name %in% input$species_select)
      }

      # Truncate long names for display
      ss$label <- ifelse(
        nchar(ss$scientific_name) > 28,
        paste0(substr(ss$scientific_name, 1, 25), "..."),
        ss$scientific_name
      )

      p <- plotly::plot_ly(ss,
        x           = ~n_passing,
        y           = ~stats::reorder(label, n_passing),
        type        = "bar",
        orientation = "h",
        name        = "Above threshold",
        marker      = list(color = "#2196F3")
      )

      below <- ss |> dplyr::mutate(n_below = n_queries - n_passing) |>
        dplyr::filter(n_below > 0L)

      if (nrow(below) > 0L) {
        p <- p |> plotly::add_trace(
          data        = below,
          x           = ~n_below,
          y           = ~stats::reorder(label, n_passing),
          type        = "bar",
          orientation = "h",
          name        = "Below threshold",
          marker      = list(color = "rgba(33,150,243,0.25)")
        )
      }

      p |>
        plotly::layout(
          barmode = "stack",
          xaxis   = list(title = "Queries"),
          yaxis   = list(title = ""),
          legend  = list(orientation = "h", y = -0.2),
          margin  = list(l = 160)
        )
    })

    # Return reactive list consumed by mod_map and DT table
    shiny::reactive({
      la <- lca_all()
      shiny::req(!is.null(la))

      filtered <- la |>
        dplyr::filter(
          best_pident >= input$min_pident,
          best_qcov   >= input$min_qcov,
          !is.na(lca_taxon),
          lca_taxon   != "unresolved",
          lca_rank    == "species"
        )

      sel <- if (length(input$species_select) > 0L) input$species_select else NULL
      if (!is.null(sel)) {
        filtered <- filtered |> dplyr::filter(lca_taxon %in% sel)
      }

      list(
        lca              = filtered,
        selected_species = unique(filtered$lca_taxon),
        thresholds       = list(min_qcov   = input$min_qcov,
                                min_pident = input$min_pident)
      )
    })
  })
}
