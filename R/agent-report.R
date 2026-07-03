#' Summarize enrichment agent results
#'
#' @param result_table Enrichment result table.
#' @param agent_name Name of the agent.
#' @param result_type Type of result being summarized.
#' @param top_n Number of top terms to include.
#'
#' @return Character summary.
#' @export
summarize_agent_results <- function(
  result_table,
  agent_name,
  result_type,
  top_n = 10
) {
  if (is.null(result_table) || nrow(result_table) == 0) {
    return(paste0(
      agent_name,
      " did not identify significant enriched ",
      result_type,
      "."
    ))
  }

  top_terms <- head(result_table$Description, top_n)

  paste0(
    agent_name, " identified enriched ", result_type, ".\n\n",
    "Top results:\n",
    paste0(seq_along(top_terms), ". ", top_terms, collapse = "\n")
  )
}

#' Generate combined KEGG and GO agent report
#'
#' @param kegg_result Output from run_kegg_agent().
#' @param go_result Output from run_go_agent().
#' @param output_file Path to output Markdown report.
#'
#' @return Report text invisibly.
#' @export
generate_two_agent_report <- function(
  kegg_result,
  go_result,
  output_file = "output/combined_report.md"
) {
  report <- paste0(
    "# Two-Agent Onco-IAN Prototype Report\n\n",
    "## Input Summary\n\n",
    "- Input genes: ", length(kegg_result$input_genes), "\n",
    "- Mapped genes: ", length(unique(kegg_result$mapped_genes$ENTREZID)), "\n\n",
    "## KEGG Agent Summary\n\n",
    kegg_result$summary,
    "\n\n",
    "## GO Agent Summary\n\n",
    go_result$summary,
    "\n\n",
    "## Combined Interpretation\n\n",
    "The KEGG agent provides pathway-level interpretation, while the GO agent provides biological process-level interpretation. ",
    "Together, these two agents provide complementary views of the DEG list and form the first working prototype of an Onco-IAN style multi-agent system.\n\n",
    "## Next Steps\n\n",
    "- Review top pathways and biological processes with the biology team.\n",
    "- Add Reactome, STRING, WikiPathways, and ChEA agents later using the same modular pattern.\n",
    "- Replace or extend rule-based summaries with LLM-generated grounded interpretations.\n"
  )

  writeLines(report, output_file)
  invisible(report)
}
