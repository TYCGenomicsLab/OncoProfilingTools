#' Run Reactome pathway enrichment agent
#'
#' Maps human gene symbols to Entrez identifiers and performs Reactome
#' pathway enrichment analysis.
#'
#' @param genes Character vector of human gene symbols.
#' @param pvalue_cutoff Numeric p-value cutoff.
#' @param p_adjust_method Multiple-testing adjustment method.
#' @param min_gene_set_size Minimum pathway size.
#' @param max_gene_set_size Maximum pathway size.
#' @param readable Logical; convert Entrez identifiers to gene symbols.
#'
#' @return A list containing Reactome results, mapped genes, and summary.
#' @export
run_reactome_agent <- function(
  genes,
  pvalue_cutoff = 0.05,
  p_adjust_method = "BH",
  min_gene_set_size = 10,
  max_gene_set_size = 500,
  readable = TRUE
) {
  if (missing(genes) || length(genes) == 0L) {
    stop("genes must be a non-empty character vector.", call. = FALSE)
  }

  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    stop("ReactomePA is required to run the Reactome agent.", call. = FALSE)
  }

  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0L) {
    stop("No valid gene symbols remained after preprocessing.", call. = FALSE)
  }

  mapped_genes <- suppressMessages(
    clusterProfiler::bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  )

  if (is.null(mapped_genes) || nrow(mapped_genes) == 0L) {
    stop("No genes could be mapped to ENTREZID.", call. = FALSE)
  }

  entrez_genes <- unique(as.character(mapped_genes$ENTREZID))

  reactome_result <- ReactomePA::enrichPathway(
    gene = entrez_genes,
    organism = "human",
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    readable = readable
  )

  reactome_table <- as.data.frame(reactome_result)

  summary <- summarize_agent_results(
    reactome_table,
    agent_name = "Reactome Agent",
    result_type = "Reactome pathways"
  )

  list(
    agent_name = "Reactome",
    input_genes = genes,
    mapped_genes = mapped_genes,
    results = reactome_table,
    enrichment_object = reactome_result,
    summary = summary
  )
}
