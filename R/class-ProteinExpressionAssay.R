# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' ProteinExpressionAssay S4 class definition
#'
#' `ProteinExpressionAssay` is a `SummarizedExperiment` subclass for storing
#' protein-coding gene expression data.
#'
#' The primary assay is named `"protein_expression"` and should contain a
#' numeric matrix-like object with features as rows and samples/models as
#' columns.
#'
#' @details
#' This class follows the Bioconductor `SummarizedExperiment` convention:
#'
#' - rows represent features, such as protein-coding genes
#' - columns represent samples, models, or cell lines
#' - `rowData()` stores feature-level metadata
#' - `colData()` stores sample/model-level metadata
#' - `metadata()` stores assay-level metadata such as units and provenance
#'
#' @return An object of class `ProteinExpressionAssay`.
#' @name ProteinExpressionAssay-class
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
methods::setClass(
  "ProteinExpressionAssay",
  contains = "SummarizedExperiment"
)

#' Validity method for ProteinExpressionAssay
#'
#' Ensures that the object contains a numeric assay named
#' `"protein_expression"` with row and column names.
#' @name internal-ProteinExpressionAssay-validity
#' @keywords internal
S4Vectors::setValidity2("ProteinExpressionAssay", function(object) {
  messages <- character()

  assay_names <- SummarizedExperiment::assayNames(object)

  if (!"protein_expression" %in% assay_names) {
    messages <- c(
      messages,
      "`ProteinExpressionAssay` must contain an assay named `protein_expression`."
    )
  }

  if ("protein_expression" %in% assay_names) {
    expression_data <- SummarizedExperiment::assay(
      object,
      "protein_expression",
      withDimnames = FALSE
    )

    if (!is.numeric(expression_data)) {
      messages <- c(
        messages,
        "`ProteinExpressionAssay` assay `protein_expression` must be numeric."
      )
    }
  }

  if (is.null(rownames(object))) {
    messages <- c(
      messages,
      "`ProteinExpressionAssay` must have row names for genes/features."
    )
  }

  if (is.null(colnames(object))) {
    messages <- c(
      messages,
      "`ProteinExpressionAssay` must have column names for samples/models."
    )
  }

  if (length(messages) > 0L) {
    messages
  } else {
    TRUE
  }
})

#' Create a ProteinExpressionAssay object
#'
#' Creates a `ProteinExpressionAssay`, a `SummarizedExperiment` subclass for
#' protein-coding gene expression data.
#'
#' @param data A numeric matrix-like object with features as rows and
#'   samples/models as columns.
#' @param rowData Optional feature-level metadata. Must have one row per feature.
#' @param colData Optional sample/model-level metadata. Must have one row per
#'   sample/model.
#' @param metadata Optional assay-level metadata list.
#' @param assay_name Name of the primary assay. Defaults to
#'   `"protein_expression"`.
#' @param unit Character string describing the expression unit, such as
#'   `"log2(TPM+1)"`.
#' @param normalized Logical scalar indicating whether the expression data are
#'   normalized.
#' @param feature_type Character string describing the feature type, such as
#'   `"protein_coding_gene"`.
#' @param source Optional character string describing the data source.
#' @param source_file Optional character string storing the source file path.
#'
#' @return A `ProteinExpressionAssay` object.
#' @name ProteinExpressionAssay
#' @export
ProteinExpressionAssay <- function(
  data,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  assay_name = "protein_expression",
  unit = "log2(TPM+1)",
  normalized = TRUE,
  feature_type = "protein_coding_gene",
  source = NA_character_,
  source_file = NA_character_
) {
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  if (!is.numeric(data)) {
    stop("`data` must be numeric.", call. = FALSE)
  }

  if (is.null(rownames(data))) {
    stop("`data` must have row names for genes/features.", call. = FALSE)
  }

  if (is.null(colnames(data))) {
    stop("`data` must have column names for samples/models.", call. = FALSE)
  }

  if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    stop("`assay_name` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit)) {
    stop("`unit` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.logical(normalized) || length(normalized) != 1L || is.na(normalized)) {
    stop("`normalized` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.character(feature_type) || length(feature_type) != 1L || is.na(feature_type) || !nzchar(feature_type)) {
    stop("`feature_type` must be a single non-empty character string.", call. = FALSE)
  }

  if (is.null(rowData)) {
    rowData <- S4Vectors::DataFrame(
      feature_id = rownames(data),
      feature_type = feature_type,
      row.names = rownames(data)
    )
  } else {
    rowData <- S4Vectors::DataFrame(rowData)

    if (nrow(rowData) != nrow(data)) {
      stop(
        "`rowData` must have the same number of rows as `data`.",
        call. = FALSE
      )
    }

    rownames(rowData) <- rownames(data)
  }

  if (is.null(colData)) {
    colData <- S4Vectors::DataFrame(
      sample_id = colnames(data),
      row.names = colnames(data)
    )
  } else {
    colData <- S4Vectors::DataFrame(colData)

    if (nrow(colData) != ncol(data)) {
      stop(
        "`colData` must have the same number of rows as columns in `data`.",
        call. = FALSE
      )
    }

    rownames(colData) <- colnames(data)
  }

  metadata <- c(
    metadata,
    list(
      unit = unit,
      normalized = normalized,
      feature_type = feature_type,
      source = source,
      source_file = source_file
    )
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(
      protein_expression = data
    ),
    rowData = rowData,
    colData = colData,
    metadata = metadata
  )

  methods::as(se, "ProteinExpressionAssay")
}
