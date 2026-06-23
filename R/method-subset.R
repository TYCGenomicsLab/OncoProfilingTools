# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Subset an OncoExperiment
#'
#' Subset an `OncoExperiment` by evaluating the expression against the
#' experiment-level `colData()`. This keeps the assay columns and metadata in
#' sync while filtering models.
#'
#' @param object An `OncoExperiment` object to subset.
#' @param subset A logical expression evaluated in the context of `colData(object)`.
#' @param select Currently unused. Included for compatibility with the base `subset()` generic.
#' @return A new `OncoExperiment` object containing only the selected models.
#' @export
methods::setMethod(
  f = "subset",
  signature = c(
    object = "OncoExperiment",
    subset = "ANY",
    select = "ANY"
  ),
  definition = function(
    object,
    subset,
    select = TRUE,
    ...
  ) {
    if (!is(object, "OncoExperiment")) {
      stop("`object` must be an OncoExperiment.", call. = FALSE)
    }

    validObject(object)

    col_data <- MultiAssayExperiment::colData(object)

    if (nrow(col_data) == 0L) {
      stop("Cannot subset because `colData(object)` is empty.", call. = FALSE)
    }

    if (ncol(col_data) == 0L) {
      stop("Cannot subset because `colData(object)` has no columns.", call. = FALSE)
    }

    col_data_df <- as.data.frame(col_data, stringsAsFactors = FALSE)

    if (missing(subset)) {
      keep <- rep(TRUE, nrow(col_data_df))
    } else {
      subset_expr <- substitute(subset)
      subset_env <- list2env(col_data_df, parent = parent.frame())

      keep <- tryCatch(
        eval(subset_expr, envir = subset_env, enclos = parent.frame()),
        error = function(e) {
          stop("Error in subset expression: ", e$message, call. = FALSE)
        }
      )
    }

    if (!is.logical(keep)) {
      stop("Result of subset expression must be a logical vector.", call. = FALSE)
    }

    if (length(keep) != nrow(col_data_df)) {
      stop(
        "Result of subset expression has incorrect length. Expected ",
        nrow(col_data_df),
        ", got ",
        length(keep),
        ".",
        call. = FALSE
      )
    }

    keep <- keep & !is.na(keep)
    keep_ids <- rownames(col_data_df)[keep]

    if (length(keep_ids) == 0L) {
      warning(
        "Subset expression matched zero models/samples. Returning an empty OncoExperiment.",
        call. = FALSE
      )
    }

    out <- object[, keep_ids, drop = FALSE]
    validObject(out)
    out
  }
)

#' Extract a column from OncoExperiment metadata
#'
#' This is a convenience accessor that forwards `$` lookups to the top-level
#' `colData()` stored on the `OncoExperiment`.
#'
#' @keywords internal
methods::setMethod(
  f = "$",
  signature = "OncoExperiment",
  definition = function(x, name) {
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
      stop("`name` must be a single non-empty character string.", call. = FALSE)
    }

    col_data <- MultiAssayExperiment::colData(x)

    if (name %in% colnames(col_data)) {
      return(col_data[[name]])
    }

    NULL
  }
)
