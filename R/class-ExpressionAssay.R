# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' ExpressionAssay S4 class definition
#'
#' `ExpressionAssay` is a `SummarizedExperiment` subclass for storing gene or
#' transcript expression data.
#'
#' The primary assay is named `"expression"` and should contain a numeric
#' matrix-like object with features as rows and samples/models as columns.
#'
#' @details
#' This class follows the Bioconductor `SummarizedExperiment` convention:
#'
#' - rows represent features, such as genes or transcripts
#' - columns represent samples, models, or cell lines
#' - `rowData()` stores feature-level metadata
#' - `colData()` stores sample/model-level metadata
#' - `metadata()` stores assay-level metadata such as units and provenance
#'
#' @return An object of class `ExpressionAssay`.
#' @name ExpressionAssay-class
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
methods::setClass(
  "ExpressionAssay",
  contains = "SummarizedExperiment"
)

#' Validity method for ExpressionAssay
#'
#' Ensures that the object contains a numeric assay named `"expression"` with
#' row and column names.
#' @name internal-ExpressionAssay-validity
#' @keywords internal
S4Vectors::setValidity2("ExpressionAssay", function(object) {
  messages <- character()

  assay_names <- SummarizedExperiment::assayNames(object)

  if (!"expression" %in% assay_names) {
    messages <- c(
      messages,
      "`ExpressionAssay` must contain an assay named `expression`."
    )
  }

  if ("expression" %in% assay_names) {
    expression_data <- SummarizedExperiment::assay(
      object,
      "expression",
      withDimnames = FALSE
    )

    if (!is.numeric(expression_data)) {
      messages <- c(
        messages,
        "`ExpressionAssay` assay `expression` must be numeric."
      )
    }
  }

  if (is.null(rownames(object))) {
    messages <- c(
      messages,
      "`ExpressionAssay` must have row names for genes/features."
    )
  }

  if (is.null(colnames(object))) {
    messages <- c(
      messages,
      "`ExpressionAssay` must have column names for samples/models."
    )
  }

  if (length(messages) > 0L) {
    messages
  } else {
    TRUE
  }
})

#' Create an ExpressionAssay object
#'
#' Creates an `ExpressionAssay`, a `SummarizedExperiment` subclass for expression
#' data.
#'
#' @param data A numeric matrix-like object with features as rows and
#'   samples/models as columns.
#' @param rowData Optional feature-level metadata. Must have one row per feature.
#' @param colData Optional sample/model-level metadata. Must have one row per
#'   sample/model.
#' @param metadata Optional assay-level metadata list.
#' @param assay_name Name of the primary assay. Defaults to `"expression"`.
#' @param unit Character string describing the expression unit, such as
#'   `"log2(TPM+1)"`.
#' @param normalized Logical scalar indicating whether the expression data are
#'   normalized.
#' @param feature_type Character string describing the feature type, such as
#'   `"gene"` or `"transcript"`.
#' @param source Optional character string describing the data source.
#' @param source_file Optional character string storing the source file path.
#'
#' @return An `ExpressionAssay` object.
#'
#' @export
ExpressionAssay <- function(
  data,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  assay_name = "expression",
  unit = "log2(TPM+1)",
  normalized = TRUE,
  feature_type = "transcript",
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
      expression = data
    ),
    rowData = rowData,
    colData = colData,
    metadata = metadata
  )

  methods::as(se, "ExpressionAssay")
}