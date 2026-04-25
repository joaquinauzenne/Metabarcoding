#' Detect and normalise BLAST result file format
#'
#' Reads the first line of a file and classifies it as:
#' \describe{
#'   \item{web}{NCBI BLAST web interface CSV — has header row containing
#'     "Per. Ident" or "Per. ident", comma-delimited, percent signs in
#'     Query Cover and Per. ident columns.}
#'   \item{cli}{Command-line blastn format 6 TSV — tab-delimited, no header,
#'     10 columns: qseqid sseqid pident qcovs evalue bitscore staxids
#'     sscinames scomnames stitle}
#' }
#'
#' @param path Path to the uploaded file (character).
#' @return A named list:
#'   \describe{
#'     \item{format}{"web" or "cli"}
#'     \item{data}{Standardised data.frame with canonical columns}
#'   }
#' @noRd
detect_and_parse_blast <- function(path) {
  first_line <- readLines(path, n = 1L, warn = FALSE)

  if (grepl("Per\\. [Ii]dent|Query [Cc]over|Scientific Name", first_line,
            ignore.case = FALSE)) {
    list(format = "web",  data = parse_blast_web(path))
  } else {
    list(format = "cli",  data = parse_blast_cli(path))
  }
}
