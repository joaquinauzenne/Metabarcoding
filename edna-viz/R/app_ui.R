#' App UI
#'
#' shinydashboard layout:
#'   Sidebar — upload module, filter module
#'   Main body — tabbed: Species Table, Map, Comparison
#'
#' @noRd
app_ui <- function(request) {
  shinydashboard::dashboardPage(
    skin = "blue",

    header = shinydashboard::dashboardHeader(
      title = "eDNA Metabarcoding Viz"
    ),

    sidebar = shinydashboard::dashboardSidebar(
      width = 300,
      shiny::div(
        class = "p-3",
        shiny::h4("Data Input", class = "sidebar-title"),
        mod_upload_ui("upload_1"),
        mod_filter_ui("filter_1")
      )
    ),

    body = shinydashboard::dashboardBody(
      # Bootstrap / custom CSS
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
          .sidebar-title { color: #ecf0f1; font-size: 14px; text-transform: uppercase;
                           letter-spacing: 1px; margin-bottom: 8px; }
          .content-wrapper { background-color: #f4f6f9; }
          .badge { font-size: 12px; padding: 4px 8px; }
          .leaflet-container { border-radius: 4px; }
        "))
      ),

      shinydashboard::tabBox(
        width = 12,

        # Tab 1: Species table
        shiny::tabPanel(
          title = shiny::tagList(shiny::icon("table"), " Species Table"),
          shiny::div(
            class = "p-2",
            shiny::uiOutput("table_summary"),
            DT::dataTableOutput("species_table")
          )
        ),

        # Tab 2: Map
        shiny::tabPanel(
          title = shiny::tagList(shiny::icon("map"), " Map"),
          mod_map_ui("map_1")
        ),

        # Tab 3: Comparison (placeholder for asv_otu_comparison output)
        shiny::tabPanel(
          title = shiny::tagList(shiny::icon("chart-bar"), " Comparison"),
          shiny::div(
            class = "p-3",
            shiny::p(
              "Upload a ", shiny::code("species_concordance.csv"),
              " produced by ", shiny::code("scripts/asv_otu_comparison.R"),
              " to visualise OTU vs ASV concordance."
            ),
            shiny::fileInput("concordance_file",
                             "Upload concordance CSV",
                             accept = ".csv"),
            DT::dataTableOutput("concordance_table")
          )
        )
      )
    )
  )
}
