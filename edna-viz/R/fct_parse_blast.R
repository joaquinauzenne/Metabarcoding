#' Parse BLAST results into a canonical data.frame
#'
#' Two parsers handle the two formats produced by this project:
#'
#' **Web format** (`parse_blast_web`):
#' NCBI BLAST web interface CSV with columns:
#'   Description, Scientific Name, Common Name, Taxid,
#'   Max Score, Total Score, Query Cover, E value, Per. ident,
#'   Acc. Len, Accession
#'
#' **CLI format** (`parse_blast_cli`):
#' blastn -outfmt 6 TSV with 10 tab-delimited columns (no header):
#'   qseqid sseqid pident qcovs evalue bitscore staxids sscinames scomnames stitle
#'
#' Both parsers return a data.frame with canonical columns:
#'   query_id, subject_id, pident, qcov, evalue, bitscore,
#'   taxid, scientific_name, common_name, description
#'
#' @noRd

# ── Canonical column set ──────────────────────────────────────────────────────
CANONICAL_COLS <- c(
  "query_id", "subject_id", "pident", "qcov",
  "evalue", "bitscore", "taxid", "scientific_name", "common_name", "description"
)

#' @param path Path to web-format CSV
#' @return Canonical data.frame
#' @noRd
parse_blast_web <- function(path) {
  raw <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal")

  # Normalise column names (handle minor header variations)
  names(raw) <- trimws(names(raw))
  col_map <- c(
    query_id        = NA_character_,   # web format has no query ID column
    subject_id      = "Accession",
    pident          = "Per. ident",
    qcov            = "Query Cover",
    evalue          = "E value",
    bitscore        = "Max Score",
    taxid           = "Taxid",
    scientific_name = "Scientific Name",
    common_name     = "Common Name",
    description     = "Description"
  )

  out <- data.frame(matrix(ncol = length(CANONICAL_COLS), nrow = nrow(raw)))
  names(out) <- CANONICAL_COLS

  for (canonical in names(col_map)) {
    src <- col_map[[canonical]]
    if (!is.na(src) && src %in% names(raw)) {
      out[[canonical]] <- raw[[src]]
    }
  }

  # query_id is not available in web format — use row index as placeholder
  out$query_id <- paste0("hit_", seq_len(nrow(out)))

  # Strip % signs and coerce to numeric
  out$pident <- as.numeric(sub("%", "", out$pident, fixed = TRUE))
  out$qcov   <- as.numeric(sub("%", "", out$qcov,   fixed = TRUE))
  out$evalue  <- suppressWarnings(as.numeric(out$evalue))
  out$bitscore <- suppressWarnings(as.numeric(out$bitscore))
  out$taxid   <- suppressWarnings(as.integer(out$taxid))

  # Drop completely empty rows
  out <- out[!is.na(out$scientific_name) & out$scientific_name != "", ]
  out
}

#' @param path Path to CLI-format TSV (no header, 10 columns)
#' @return Canonical data.frame
#' @noRd
parse_blast_cli <- function(path) {
  cli_cols <- c("query_id", "subject_id", "pident", "qcov",
                "evalue", "bitscore", "taxid_raw", "scientific_name",
                "common_name", "description")

  raw <- readr::read_tsv(path, col_names = cli_cols,
                          show_col_types = FALSE,
                          col_types = readr::cols(.default = readr::col_character()))

  # Handle semicolon-delimited taxids/scientific names (multiple taxa per hit)
  # Take the first entry in each semicolon-separated list
  raw$taxid           <- suppressWarnings(
    as.integer(sub(";.*", "", raw$taxid_raw))
  )
  raw$scientific_name <- trimws(sub(";.*", "", raw$scientific_name))
  raw$common_name     <- trimws(raw$common_name)

  out <- raw[, CANONICAL_COLS, drop = FALSE]
  out$pident   <- suppressWarnings(as.numeric(out$pident))
  out$qcov     <- suppressWarnings(as.numeric(out$qcov))
  out$evalue   <- suppressWarnings(as.numeric(out$evalue))
  out$bitscore <- suppressWarnings(as.numeric(out$bitscore))

  # Drop empty rows
  out <- out[!is.na(out$scientific_name) & out$scientific_name != "", ]
  out
}

#' Validate canonical BLAST data.frame
#'
#' Returns a character vector of warning messages (empty if valid).
#' @noRd
validate_blast_data <- function(df) {
  msgs <- character(0)
  missing <- setdiff(CANONICAL_COLS, names(df))
  if (length(missing)) {
    msgs <- c(msgs, paste("Missing columns:", paste(missing, collapse = ", ")))
  }
  if (nrow(df) == 0) {
    msgs <- c(msgs, "No rows after parsing — check file format.")
  }
  if (all(is.na(df$pident))) {
    msgs <- c(msgs, "pident column is all NA — percent identity not parsed correctly.")
  }
  msgs
}
