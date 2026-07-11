#' Run GSVA pathway activity agent
#'
#' @param expression_matrix Numeric matrix with genes as rows and samples
#'   as columns. Row names must contain gene identifiers.
#' @param gene_sets Named list of gene sets. Each element must contain
#'   gene identifiers matching the row names of expression_matrix.
#' @param kcdf Kernel used by GSVA. Use "Gaussian" for continuous
#'   normalized values such as log-CPM or log-TPM.
#' @param min_size Minimum gene-set size.
#' @param max_size Maximum gene-set size.
#'
#' @return A list containing GSVA scores and analysis metadata.
#' @export
run_gsva_agent <- function(
  expression_matrix,
  gene_sets,
  kcdf = "Gaussian",
  min_size = 5,
  max_size = 500
) {
  if (missing(expression_matrix) || is.null(expression_matrix)) {
    stop("expression_matrix is required.", call. = FALSE)
  }

  if (missing(gene_sets) || length(gene_sets) == 0) {
    stop("gene_sets must be a non-empty named list.", call. = FALSE)
  }

  expression_matrix <- as.matrix(expression_matrix)

  if (!is.numeric(expression_matrix)) {
    stop("expression_matrix must contain numeric values.", call. = FALSE)
  }

  if (is.null(rownames(expression_matrix))) {
    stop(
      "expression_matrix must have gene identifiers as row names.",
      call. = FALSE
    )
  }

  if (is.null(colnames(expression_matrix))) {
    stop(
      "expression_matrix must have sample identifiers as column names.",
      call. = FALSE
    )
  }

  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("gene_sets must be a named list.", call. = FALSE)
  }

  expression_matrix <- expression_matrix[
    !duplicated(rownames(expression_matrix)),
    ,
    drop = FALSE
  ]

  gene_sets <- lapply(gene_sets, function(x) {
    unique(as.character(x))
  })

  gsva_param <- GSVA::gsvaParam(
    expression_matrix,
    gene_sets,
    kcdf = kcdf,
    minSize = min_size,
    maxSize = max_size
  )

  gsva_scores <- GSVA::gsva(
    gsva_param,
    verbose = FALSE
  )

  gsva_scores <- as.matrix(gsva_scores)

  summary <- paste0(
    "GSVA Agent calculated pathway activity scores for ",
    nrow(gsva_scores),
    " gene sets across ",
    ncol(gsva_scores),
    " samples."
  )

  list(
    agent_name = "GSVA",
    input_dimensions = dim(expression_matrix),
    gene_set_count = length(gene_sets),
    results = gsva_scores,
    summary = summary
  )
}