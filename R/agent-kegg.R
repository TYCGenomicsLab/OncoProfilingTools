#' Run KEGG enrichment agent
#'
#' @param genes Character vector of gene symbols.
#' @param organism KEGG organism code. Default is "hsa" for human.
#' @param pvalue_cutoff P-value cutoff.
#' @param max_attempts Maximum attempts for transient KEGG REST failures.
#' @param retry_delay_seconds Initial retry delay; subsequent delays use exponential backoff.
#' @param request_timeout_seconds Per-attempt network timeout in seconds. The
#'   default allows for the comparatively large KEGG pathway-link response.
#' @param enrich_fn Enrichment function. Defaults to `clusterProfiler::enrichKEGG`;
#'   injectable for deterministic tests.
#' @param sleep_fn Delay function. Defaults to `Sys.sleep`; injectable for tests.
#'
#' @return A list containing KEGG results, mapped genes, and summary.
#' @export
run_kegg_agent <- function(
  genes,
  organism = "hsa",
  pvalue_cutoff = 0.05,
  max_attempts = 3L,
  retry_delay_seconds = 2,
  request_timeout_seconds = 180,
  enrich_fn = NULL,
  sleep_fn = Sys.sleep
) {
  if (missing(genes) || length(genes) == 0) {
    stop("genes must be a non-empty character vector.", call. = FALSE)
  }

  genes <- unique(as.character(genes))

  mapped_genes <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )

  entrez_genes <- unique(mapped_genes$ENTREZID)

  if (length(entrez_genes) == 0) {
    stop("No genes could be mapped to ENTREZID.", call. = FALSE)
  }

  kegg_result <- run_kegg_enrichment_with_retry(
    gene = entrez_genes,
    organism = organism,
    pvalue_cutoff = pvalue_cutoff,
    max_attempts = max_attempts,
    retry_delay_seconds = retry_delay_seconds,
    request_timeout_seconds = request_timeout_seconds,
    enrich_fn = enrich_fn,
    sleep_fn = sleep_fn
  )

  kegg_table <- as.data.frame(kegg_result)

  summary <- summarize_agent_results(
    kegg_table,
    agent_name = "KEGG Agent",
    result_type = "KEGG pathways"
  )

  list(
    agent_name = "KEGG",
    input_genes = genes,
    mapped_genes = mapped_genes,
    results = kegg_table,
    summary = summary
  )
}

is_transient_kegg_error <- function(error) {
  grepl(
    paste(
      "cannot read from connection",
      "timed? ?out|timeout",
      "connection (reset|refused|closed)",
      "could not resolve|couldn't connect|failed to connect",
      "temporary failure|server returned nothing",
      "http.*(429|500|502|503|504)",
      sep = "|"
    ),
    conditionMessage(error),
    ignore.case = TRUE,
    perl = TRUE
  )
}

run_kegg_enrichment_with_retry <- function(
  gene,
  organism = "hsa",
  pvalue_cutoff = 0.05,
  max_attempts = 3L,
  retry_delay_seconds = 2,
  request_timeout_seconds = 180,
  enrich_fn = NULL,
  sleep_fn = Sys.sleep
) {
  max_attempts <- suppressWarnings(as.integer(max_attempts))
  if (!is.finite(max_attempts) || max_attempts < 1L) max_attempts <- 1L
  retry_delay_seconds <- suppressWarnings(as.numeric(retry_delay_seconds))
  if (!is.finite(retry_delay_seconds) || retry_delay_seconds < 0) retry_delay_seconds <- 0
  request_timeout_seconds <- suppressWarnings(as.numeric(request_timeout_seconds))
  if (!is.finite(request_timeout_seconds) || request_timeout_seconds < 1) request_timeout_seconds <- 180
  if (is.null(enrich_fn)) enrich_fn <- clusterProfiler::enrichKEGG

  previous_timeout <- getOption("timeout", 60)
  options(timeout = max(as.numeric(previous_timeout), request_timeout_seconds))
  on.exit(options(timeout = previous_timeout), add = TRUE)

  last_error <- NULL
  attempts_used <- 0L
  for (attempt in seq_len(max_attempts)) {
    attempts_used <- attempt
    value <- tryCatch(
      enrich_fn(
        gene = gene,
        organism = organism,
        pvalueCutoff = pvalue_cutoff
      ),
      error = function(error) error
    )
    if (!inherits(value, "condition")) return(value)

    last_error <- value
    transient <- is_transient_kegg_error(value)
    if (!transient || attempt >= max_attempts) break

    delay <- retry_delay_seconds * (2 ^ (attempt - 1L))
    message(
      "KEGG REST request attempt ", attempt, " of ", max_attempts,
      " failed transiently (", conditionMessage(value), "). Retrying in ",
      format(delay, trim = TRUE), " seconds."
    )
    if (delay > 0) sleep_fn(delay)
  }

  stop(
    paste0(
      "KEGG REST enrichment failed after ", attempts_used,
      if (attempts_used == 1L) " attempt: " else " attempts: ",
      conditionMessage(last_error),
      ". This is a KEGG service/network failure; mapped genes remain valid. Retry when rest.kegg.jp is reachable."
    ),
    call. = FALSE
  )
}
