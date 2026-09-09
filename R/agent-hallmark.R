#' Run Hallmark gene-set enrichment agent
#'
#' Performs over-representation analysis using the curated Hallmark
#' gene-set collection from MSigDB.
#'
#' @param genes Character vector of human gene symbols.
#' @param pvalue_cutoff P-value cutoff. Default is 0.05.
#' @param qvalue_cutoff Q-value cutoff. Default is 0.20.
#'
#' @return A list containing Hallmark enrichment results,
#' input genes, mapped genes, and a text summary.
#' @export
run_hallmark_agent <- function(
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

  hallmark_sets <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "H"
  )

  term2gene <- unique(
    hallmark_sets[, c("gs_name", "gene_symbol")]
  )

  colnames(term2gene) <- c("term", "gene")

  mapped_genes <- intersect(
    genes,
    unique(term2gene$gene)
  )

  if (length(mapped_genes) == 0) {
    stop(
      "None of the supplied genes were found in the Hallmark gene sets.",
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
    agent_name = "Hallmark Agent",
    result_type = "Hallmark gene sets"
  )

  list(
    agent_name = "Hallmark",
    input_genes = genes,
    mapped_genes = mapped_genes,
    results = result_table,
    summary = summary
  )
}
