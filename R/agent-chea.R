#' Run ChEA transcription factor enrichment agent
#'
#' @param genes Character vector of human gene symbols.
#' @param database Enrichr ChEA database name.
#'
#' @return A list containing ChEA enrichment results and summary.
#' @export
run_chea_agent <- function(
  genes,
  database = "ChEA_2022"
) {
  if (missing(genes) || length(genes) == 0) {
    stop(
      "genes must be a non-empty character vector.",
      call. = FALSE
    )
  }

  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop(
      "Package 'enrichR' is required.",
      call. = FALSE
    )
  }

  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & genes != ""]

  if (length(genes) == 0) {
    stop(
      "No valid genes remained after cleaning.",
      call. = FALSE
    )
  }

  enrichment <- enrichR::enrichr(
    genes,
    database
  )

  if (!database %in% names(enrichment)) {
    stop(
      paste("No results returned for", database),
      call. = FALSE
    )
  }

  chea_table <- enrichment[[database]]

  summary <- if (nrow(chea_table) == 0) {
    paste(
      "ChEA Agent found no enriched transcription factors for",
      length(genes),
      "input genes."
    )
  } else {
    paste(
      "ChEA Agent identified",
      nrow(chea_table),
      "transcription-factor terms from",
      length(genes),
      "input genes."
    )
  }

  list(
    agent_name = "ChEA",
    database = database,
    input_genes = genes,
    results = chea_table,
    summary = summary
  )
}