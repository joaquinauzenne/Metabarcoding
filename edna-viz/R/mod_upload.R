#' Upload Module — File Input & Format Detection
#'
#' Accepts BLAST results in two formats:
#'   - Web CSV: NCBI BLAST web interface output
#'   - CLI TSV: blastn -outfmt 6 (10-column, produced by pipeline.sh)
#'
#' Also provides a "Load demo" button to load the bundled example file.
#'
#' @param id Shiny module ID.
#' @return (server) Reactive returning canonical BLAST data.frame, or NULL.

# ── UI ────────────────────────────────────────────────────────────────────────
mod_upload_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fileInput(
      ns("file"),
      label       = "Upload BLAST results",
      accept      = c(".csv", ".tsv", ".txt"),
      buttonLabel = "Browse...",
      placeholder = "Web CSV or CLI TSV"
    ),
    shiny::actionButton(
      ns("load_demo"),
      label = "Load demo (Alouatta example)",
      icon  = shiny::icon("flask"),
      class = "btn-sm btn-outline-secondary"
    ),
    shiny::br(),
    shiny::uiOutput(ns("format_badge")),
    shiny::uiOutput(ns("parse_warnings"))
  )
}

# ── Server ────────────────────────────────────────────────────────────────────
mod_upload_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Internal reactive: raw parse result (list with $format and $data)
    parsed <- shiny::reactiveVal(NULL)

    # Upload from file input
    shiny::observeEvent(input$file, {
      req(input$file)
      result <- tryCatch(
        detect_and_parse_blast(input$file$datapath),
        error = function(e) list(error = conditionMessage(e))
      )
      parsed(result)
    })

    # Load bundled demo file
    shiny::observeEvent(input$load_demo, {
      demo_path <- system.file(
        "app/app_data/demo_blast.csv",
        package = "ednaviz"
      )
      # Fall back to examples/ directory relative to package root
      if (!file.exists(demo_path)) {
        demo_path <- system.file(
          "../../examples/2ER9WHC1013-Alignment-Descriptions.csv",
          package = "ednaviz"
        )
      }
      if (!file.exists(demo_path)) {
        shiny::showNotification("Demo file not found.", type = "error")
        return()
      }
      result <- tryCatch(
        detect_and_parse_blast(demo_path),
        error = function(e) list(error = conditionMessage(e))
      )
      parsed(result)
      shiny::showNotification("Demo data loaded.", type = "message", duration = 3)
    })

    # Format badge
    output$format_badge <- shiny::renderUI({
      p <- parsed()
      if (is.null(p)) return(NULL)
      if (!is.null(p$error)) {
        shiny::tags$span(class = "badge bg-danger", "Parse error")
      } else {
        label <- if (p$format == "web") "Web CSV" else "CLI TSV"
        n     <- nrow(p$data)
        shiny::tags$span(
          class = "badge bg-success",
          sprintf("%s — %d rows", label, n)
        )
      }
    })

    # Parse warnings / errors
    output$parse_warnings <- shiny::renderUI({
      p <- parsed()
      if (is.null(p)) return(NULL)
      if (!is.null(p$error)) {
        shiny::tags$div(
          class = "alert alert-danger mt-2",
          shiny::tags$strong("Error: "), p$error
        )
      } else {
        warnings <- validate_blast_data(p$data)
        if (length(warnings) == 0L) return(NULL)
        shiny::tags$div(
          class = "alert alert-warning mt-2",
          shiny::tags$ul(lapply(warnings, shiny::tags$li))
        )
      }
    })

    # Return reactive: canonical data.frame or NULL
    shiny::reactive({
      p <- parsed()
      if (is.null(p) || !is.null(p$error)) return(NULL)
      p$data
    })
  })
}
