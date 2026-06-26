# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

methods::setClass(
  "MutationAssay",
  contains = "SummarizedExperiment"
)

S4Vectors::setValidity2("MutationAssay", function(object) {
  messages <- character()

  assay_names <- SummarizedExperiment::assayNames(object)

  if (!"mutation" %in% assay_names) {
    messages <- c(
      messages,
      "`MutationAssay` must contain an assay named `mutation`."
    )
  }

  if ("mutation" %in% assay_names) {
    mutation_data <- SummarizedExperiment::assay(
      object,
      "mutation",
      withDimnames = FALSE
    )

    if (!is.numeric(mutation_data)) {
      messages <- c(
        messages,
        "`MutationAssay` assay `mutation` must be numeric."
      )
    }
  }

  if (is.null(rownames(object))) {
    messages <- c(
      messages,
      "`MutationAssay` must have row names for mutation calls."
    )
  }

  if (is.null(colnames(object))) {
    messages <- c(
      messages,
      "`MutationAssay` must have column names for samples/models."
    )
  }

  if (length(messages) > 0L) {
    messages
  } else {
    TRUE
  }
})

MutationAssay <- function(
  data,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  assay_name = "mutation"
) {
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  if (!is.numeric(data)) {
    stop("`data` must be numeric.", call. = FALSE)
  }

  if (is.null(rownames(data))) {
    stop("`data` must have row names for mutation calls.", call. = FALSE)
  }

  if (is.null(colnames(data))) {
    stop("`data` must have column names for samples/models.", call. = FALSE)
  }

  if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    stop("`assay_name` must be a single non-empty character string.", call. = FALSE)
  }

  if (is.null(rowData)) {
    rowData <- S4Vectors::DataFrame(
      mutation_call_id = rownames(data),
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
      assay_name = assay_name
    )
  )

  assays <- S4Vectors::SimpleList()
  assays[[assay_name]] <- data

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = rowData,
    colData = colData,
    metadata = metadata
  )

  methods::as(se, "MutationAssay")
}
