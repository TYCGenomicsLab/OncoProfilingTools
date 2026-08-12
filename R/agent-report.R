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
    return(paste0(agent_name, " did not identify significant enriched ", result_type, "."))
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
  input_n <- length(kegg_result$input_genes)
  mapped_n <- length(unique(kegg_result$mapped_genes$ENTREZID))

  kegg_n <- ifelse(is.null(kegg_result$results), 0, nrow(kegg_result$results))
  go_n <- ifelse(is.null(go_result$results), 0, nrow(go_result$results))

  report <- paste0(
    "# Two-Agent Onco-IAN Prototype Report\n\n",
    "## 1. Prototype Status\n\n",
    "This report was generated from the current two-agent prototype developed this week. ",
    "The system currently supports two enrichment agents: a KEGG pathway enrichment agent and a GO Biological Process enrichment agent.\n\n",
    "## 2. Input Summary\n\n",
    "- Input genes analyzed: ", input_n, "\n",
    "- Successfully mapped genes: ", mapped_n, "\n",
    "- Mapping rate: ", round((mapped_n / input_n) * 100, 2), "%\n",
    "- KEGG enriched terms returned: ", kegg_n, "\n",
    "- GO enriched terms returned: ", go_n, "\n\n",
    "## 3. KEGG Agent Output\n\n",
    kegg_result$summary,
    "\n\n",
    "## 4. GO Biological Process Agent Output\n\n",
    go_result$summary,
    "\n\n",
    "## 5. Combined Interpretation\n\n",
    "The KEGG agent provides pathway-level interpretation, while the GO agent provides biological process-level interpretation. ",
    "Together, these agents provide complementary biological views of the input DEG list. ",
    "This confirms that the initial Onco-IAN style workflow can take a gene list, map gene symbols to Entrez identifiers, run enrichment analysis, and generate a structured report.\n\n",
    "## 6. Current Limitations\n\n",
    "- The current prototype uses rule-based summaries rather than full LLM-generated interpretation.\n",
    "- The workflow currently includes KEGG and GO only.\n",
    "- Biological interpretation still needs review by the research team.\n",
    "- Dataset-specific assumptions should be validated before paper-level conclusions.\n\n",
    "## 7. Next Steps\n\n",
    "- Run the same workflow on the CRC/CMS4 dataset provided by the data team.\n",
    "- Improve report generation with more detailed interpretation sections.\n",
    "- Add additional agents such as Reactome, STRING, WikiPathways, and ChEA.\n",
    "- Prepare screenshots and results for weekly lab progress presentation.\n"
  )

  writeLines(report, output_file)
  invisible(report)
}
