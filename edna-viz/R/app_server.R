#' App Server
#' @noRd
app_server <- function(input, output, session) {

  # Module: file upload → canonical BLAST data.frame (reactive)
  blast_data <- mod_upload_server("upload_1")

  # Module: threshold filtering + LCA → filtered species list (reactive list)
  filtered_data <- mod_filter_server("filter_1", blast_data = blast_data)

  # Module: Leaflet map — consumes filtered species data
  mod_map_server("map_1", species_data = filtered_data)

  # ── Species Table (main panel) ───────────────────────────────────────────────
  output$table_summary <- shiny::renderUI({
    fd <- filtered_data()
    if (is.null(fd) || nrow(fd$lca) == 0L) {
      return(shiny::p("No species passing current thresholds.",
                      class = "text-muted"))
    }
    n_sp <- length(unique(fd$lca$lca_taxon))
    shiny::p(sprintf("%d species detected above thresholds.", n_sp),
             class = "text-success mb-2")
  })

  output$species_table <- DT::renderDataTable({
    fd <- filtered_data()
    shiny::req(!is.null(fd), nrow(fd$lca) > 0L)

    display <- fd$lca |>
      dplyr::select(
        `Scientific Name` = lca_taxon,
        `Rank`            = lca_rank,
        `Confidence`      = confidence,
        `Best % ID`       = best_pident,
        `Best Q Cov`      = best_qcov,
        `Queries`         = n_species,
        `All Hits`        = all_species
      ) |>
      dplyr::arrange(dplyr::desc(`Best % ID`))

    DT::datatable(
      display,
      options = list(
        pageLength  = 20,
        scrollX     = TRUE,
        dom         = "Bfrtip",
        buttons     = c("csv", "excel")
      ),
      rownames   = FALSE,
      extensions = "Buttons",
      class      = "compact stripe hover"
    ) |>
      DT::formatRound(c("Best % ID", "Best Q Cov"), digits = 2) |>
      DT::formatStyle(
        "Confidence",
        backgroundColor = DT::styleEqual(
          c("high", "medium", "low"),
          c("#d4edda", "#fff3cd", "#f8d7da")
        )
      )
  })

  # ── Comparison tab: species concordance CSV ───────────────────────────────
  concordance_data <- shiny::reactive({
    shiny::req(input$concordance_file)
    tryCatch(
      readr::read_csv(input$concordance_file$datapath, show_col_types = FALSE),
      error = function(e) NULL
    )
  })

  output$concordance_table <- DT::renderDataTable({
    cd <- concordance_data()
    shiny::req(!is.null(cd))
    DT::datatable(cd, rownames = FALSE,
                  options = list(scrollX = TRUE, pageLength = 25))
  })
}
