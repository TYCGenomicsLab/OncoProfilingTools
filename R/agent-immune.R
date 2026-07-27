#' Run immune deconvolution agent
#'
#' Estimates immune-cell abundance from a bulk gene-expression matrix.
#'
#' @param expression_matrix Numeric matrix or data frame with genes as rows
#'   and samples as columns. Row names must contain human gene symbols.
#' @param method Immune deconvolution method, such as "quantiseq", "epic",
#'   "mcp_counter", "xcell", or another method supported by immunedeconv.
#' @param tumor Logical indicating whether samples are tumor samples.
#' @param arrays Logical indicating whether input came from microarrays.
#'
#' @return A list containing immune-cell estimates and analysis metadata.
#' @export
run_immune_agent <- function(
  expression_matrix,
  method = "quantiseq",
  tumor = TRUE,
  arrays = FALSE
) {
  if (missing(expression_matrix) || is.null(expression_matrix)) {
    stop("expression_matrix is required.", call. = FALSE)
  }

  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    stop(
      paste0(
        "The immunedeconv package is required. Install it using ",
        "remotes::install_github('omnideconv/immunedeconv')."
      ),
      call. = FALSE
    )
  }

  expression_matrix <- as.matrix(expression_matrix)

  if (!is.numeric(expression_matrix)) {
    stop("expression_matrix must contain numeric values.", call. = FALSE)
  }

  if (is.null(rownames(expression_matrix))) {
    stop(
      "expression_matrix must have gene symbols as row names.",
      call. = FALSE
    )
  }

  if (is.null(colnames(expression_matrix))) {
    stop(
      "expression_matrix must have sample identifiers as column names.",
      call. = FALSE
    )
  }

  if (nrow(expression_matrix) == 0L || ncol(expression_matrix) == 0L) {
    stop(
      "expression_matrix must contain at least one gene and one sample.",
      call. = FALSE
    )
  }

  gene_names <- trimws(rownames(expression_matrix))

  keep_rows <- !is.na(gene_names) &
    nzchar(gene_names) &
    !duplicated(gene_names)

  expression_matrix <- expression_matrix[
    keep_rows, ,
    drop = FALSE
  ]

  rownames(expression_matrix) <- gene_names[keep_rows]

  expression_matrix[!is.finite(expression_matrix)] <- NA_real_

  if (anyNA(expression_matrix)) {
    warning(
      "Missing or non-finite expression values were replaced with zero.",
      call. = FALSE
    )
    expression_matrix[is.na(expression_matrix)] <- 0
  }

  deconvolute_arguments <- list(
    gene_expression = expression_matrix,
    method = method
  )

  supported_arguments <- names(formals(immunedeconv::deconvolute))

  if ("tumor" %in% supported_arguments) {
    deconvolute_arguments$tumor <- tumor
  }

  if ("arrays" %in% supported_arguments) {
    deconvolute_arguments$arrays <- arrays
  }

  deconvolution_result <- do.call(
    immunedeconv::deconvolute,
    deconvolute_arguments
  )

  result_table <- as.data.frame(
    deconvolution_result,
    stringsAsFactors = FALSE
  )

  summary <- paste0(
    "Immune Deconvolution Agent used the ",
    method,
    " method to estimate ",
    nrow(result_table),
    " immune or stromal cell populations across ",
    ncol(expression_matrix),
    " samples."
  )

  list(
    agent_name = "Immune Deconvolution",
    method = method,
    input_dimensions = dim(expression_matrix),
    input_genes = rownames(expression_matrix),
    results = result_table,
    summary = summary
  )
}
