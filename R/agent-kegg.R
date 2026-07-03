#' Run KEGG enrichment agent
#'
#' @param genes Character vector of gene symbols.
#' @param organism KEGG organism code. Default is "hsa" for human.
#' @param pvalue_cutoff P-value cutoff.
#'
#' @return A list containing KEGG results, mapped genes, and summary.
#' @export
run_kegg_agent <- function(
  genes,
  organism = "hsa",
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

  kegg_result <- clusterProfiler::enrichKEGG(
    gene = entrez_genes,
    organism = organism,
    pvalueCutoff = pvalue_cutoff
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
