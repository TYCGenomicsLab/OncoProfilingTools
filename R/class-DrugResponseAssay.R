# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' DrugResponseAssay S4 class definition
#'
#' `DrugResponseAssay` is a `SummarizedExperiment` subclass for storing drug
#' response data.
#'
#' The primary assay is named `"response"` and should contain a numeric
#' matrix-like object with treatments/drugs as rows and samples/models as
#' columns.
#'
#' @details
#' This class follows the Bioconductor `SummarizedExperiment` convention:
#'
#' - rows represent treatments, drugs, or compounds
#' - columns represent samples, models, or cell lines
#' - `rowData()` stores treatment/compound-level metadata
#' - `colData()` stores sample/model-level metadata
#' - `metadata()` stores assay-level metadata such as units and provenance
#'
#' Missing values are allowed because drug response matrices are often partially
#' observed.
#'
#' @return An object of class `DrugResponseAssay`.
#'
#' @name DrugResponseAssay-class
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
methods::setClass(
  "DrugResponseAssay",
  contains = "SummarizedExperiment"
)

#' Validity method for DrugResponseAssay
#'
#' Ensures that the object contains a numeric assay named `"response"` with row
#' and column names.
#'
#' @keywords internal
S4Vectors::setValidity2("DrugResponseAssay", function(object) {
  messages <- character()

  assay_names <- SummarizedExperiment::assayNames(object)

  if (!"response" %in% assay_names) {
    messages <- c(
      messages,
      "`DrugResponseAssay` must contain an assay named `response`."
    )
  }

  if ("response" %in% assay_names) {
    response_data <- SummarizedExperiment::assay(
      object,
      "response",
      withDimnames = FALSE
    )

    if (!is.numeric(response_data)) {
      messages <- c(
        messages,
        "`DrugResponseAssay` assay `response` must be numeric."
      )
    }
  }

  if (is.null(rownames(object))) {
    messages <- c(
      messages,
      "`DrugResponseAssay` must have row names for treatments/drugs."
    )
  }

  if (is.null(colnames(object))) {
    messages <- c(
      messages,
      "`DrugResponseAssay` must have column names for samples/models."
    )
  }

  if (length(messages) > 0L) {
    messages
  } else {
    TRUE
  }
})

#' Create a DrugResponseAssay object
#'
#' Creates a `DrugResponseAssay`, a `SummarizedExperiment` subclass for drug
#' response data.
#'
#' @param data A numeric matrix-like object with treatments/drugs as rows and
#'   samples/models as columns.
#' @param rowData Optional treatment/compound-level metadata. Must have one row
#'   per treatment/drug.
#' @param colData Optional sample/model-level metadata. Must have one row per
#'   sample/model.
#' @param metadata Optional assay-level metadata list.
#' @param assay_name Name of the primary assay. Defaults to `"response"`.
#' @param unit Character string describing the response unit.
#' @param source Optional character string describing the data source.
#' @param source_file Optional character string storing the source file path.
#'
#' @return A `DrugResponseAssay` object.
#'
#' @export
DrugResponseAssay <- function(
  data,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  assay_name = "response",
  unit = "logfoldchange",
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
    stop("`data` must have row names for treatments/drugs.", call. = FALSE)
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

  if (is.null(rowData)) {
    rowData <- S4Vectors::DataFrame(
      treatment_id = rownames(data),
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
      model_id = colnames(data),
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
      source = source,
      source_file = source_file,
      n_missing_values = sum(is.na(data)),
      missing_fraction = mean(is.na(data))
    )
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(
      response = data
    ),
    rowData = rowData,
    colData = colData,
    metadata = metadata
  )

  methods::as(se, "DrugResponseAssay")
}
