#' Run WikiPathways enrichment agent
#'
#' Performs over-representation analysis using the WikiPathways
#' gene-set collection from MSigDB.
#'
#' @param genes Character vector of human gene symbols.
#' @param pvalue_cutoff P-value cutoff. Default is 0.05.
#' @param qvalue_cutoff Q-value cutoff. Default is 0.20.
#'
#' @return A list containing WikiPathways enrichment results,
#' input genes, mapped genes, and a text summary.
#' @export
run_wikipathways_agent <- function(
  genes,
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.20
) {
  if (missing(genes) || length(genes) == 0) {
    stop("genes must be a non-empty character vector.", call. = FALSE)
  }

  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0) {
    stop("No valid gene symbols were supplied.", call. = FALSE)
  }

  wikipathways_sets <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:WIKIPATHWAYS"
  )

  term2gene <- unique(
    wikipathways_sets[, c("gs_name", "gene_symbol")]
  )

  colnames(term2gene) <- c("term", "gene")

  mapped_genes <- intersect(
    genes,
    unique(term2gene$gene)
  )

  if (length(mapped_genes) == 0) {
    stop(
      "None of the supplied genes were found in the WikiPathways gene sets.",
      call. = FALSE
    )
  }

  enrichment_result <- clusterProfiler::enricher(
    gene = mapped_genes,
    TERM2GENE = term2gene,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
  )

  result_table <- as.data.frame(enrichment_result)

  summary <- summarize_agent_results(
    result_table,
    agent_name = "WikiPathways Agent",
    result_type = "WikiPathways terms"
  )

  list(
    agent_name = "WikiPathways",
    input_genes = genes,
    mapped_genes = mapped_genes,
    results = result_table,
    summary = summary
  )
}
