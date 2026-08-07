#' Run ChEA transcription factor enrichment agent
#'
#' @param genes Character vector of human gene symbols.
#' @param database Enrichr ChEA database name.
#' @param max_attempts Number of fresh Enrichr requests attempted when the
#'   service returns an expired or malformed response.
#' @param retry_delay_seconds Delay between retry attempts.
#' @param request_fn Optional request function used for testing. It must accept
#'   `(genes, database)` and return the same structure as [enrichR::enrichr()].
#' @param sleep_fn Optional sleep function used for testing.
#'
#' @return A list containing ChEA enrichment results and summary.
#' @export
run_chea_agent <- function(
  genes,
  database = "ChEA_2022",
  max_attempts = 3L,
  retry_delay_seconds = 1,
  request_fn = NULL,
  sleep_fn = Sys.sleep
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

  max_attempts <- max(1L, as.integer(max_attempts))
  if (is.null(request_fn)) request_fn <- enrichR::enrichr

  required_columns <- c("Term", "Adjusted.P.value", "Combined.Score")
  last_problem <- paste("No results returned for", database)
  chea_table <- NULL

  for (attempt in seq_len(max_attempts)) {
    enrichment <- tryCatch(
      request_fn(genes, database),
      error = function(error) {
        last_problem <<- conditionMessage(error)
        NULL
      }
    )

    candidate <- if (!is.null(enrichment) && database %in% names(enrichment)) {
      enrichment[[database]]
    } else {
      NULL
    }

    valid_schema <- is.data.frame(candidate) &&
      all(required_columns %in% names(candidate))
    expired_payload <- is.data.frame(candidate) &&
      any(grepl("expired", unlist(candidate, use.names = FALSE), ignore.case = TRUE))

    if (valid_schema && !expired_payload) {
      chea_table <- candidate
      break
    }

    if (!is.null(candidate)) {
      last_problem <- paste(
        "the service returned an expired or malformed response instead of",
        "a ChEA result table"
      )
    }

    if (attempt < max_attempts && retry_delay_seconds > 0) {
      sleep_fn(retry_delay_seconds)
    }
  }

  if (is.null(chea_table)) {
    stop(
      paste0(
        "ChEA could not obtain a valid Enrichr response after ",
        max_attempts, " attempts: ", last_problem,
        ". This module was marked failed; expired API data was not saved."
      ),
      call. = FALSE
    )
  }

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
