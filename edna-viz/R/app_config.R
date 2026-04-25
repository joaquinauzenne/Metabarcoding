#' Access the Golem config
#'
#' @param value Value to retrieve from the config. If NULL, returns the full config.
#' @param config The configuration to use. Default uses the GOLEM_CONFIG_ACTIVE
#'   environment variable, falling back to "default".
#' @param use_parent If TRUE, fallback to the parent config.
#' @noRd
get_golem_config <- function(value = NULL,
                              config = Sys.getenv(
                                "GOLEM_CONFIG_ACTIVE",
                                Sys.getenv("R_CONFIG_ACTIVE", "default")
                              ),
                              use_parent = TRUE) {
  config::get(
    value      = value,
    config     = config,
    file       = system.file("golem-config.yml", package = "ednaviz"),
    use_parent = use_parent
  )
}

# Global options ───────────────────────────────────────────────────────────────
options(
  shiny.maxRequestSize = 50 * 1024^2,   # 50 MB upload limit
  shiny.useragg        = TRUE            # use ragg for rendering if available
)
