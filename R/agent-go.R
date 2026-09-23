#' Run GO Biological Process enrichment agent
#'
#' @param genes Character vector of gene symbols.
#' @param ontology GO ontology. Default is "BP".
#' @param pvalue_cutoff P-value cutoff.
#'
#' @return A list containing GO results, mapped genes, and summary.
#' @export
run_go_agent <- function(
  genes,
  ontology = "BP",
  pvalue_cutoff = 0.05
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

  go_result <- clusterProfiler::enrichGO(
    gene = entrez_genes,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = ontology,
    pvalueCutoff = pvalue_cutoff,
    readable = TRUE
  )

  go_table <- as.data.frame(go_result)

  summary <- summarize_agent_results(
    go_table,
    agent_name = "GO Agent",
    result_type = "GO Biological Process terms"
  )

  list(
    agent_name = "GO",
    input_genes = genes,
    mapped_genes = mapped_genes,
    results = go_table,
    summary = summary
  )
}
