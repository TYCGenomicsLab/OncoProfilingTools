# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Subset an OncoExperiment
#'
#' Allow the user the subset an OncoExperiment by a given condition,
#' similar to the subset function in Seurat. The condition will be evaluated
#' and the heavy lifitng of subsetting the `MultiAssayExperiment` within an `OncoExperiment`
#' will be takne care of.
#'
#' @param object An `OncoExperiment` object to be subset.
#' @param subset A logical expression indicating which rows to keep.
#' @param select A logical expression indicating which columns to keep.
#' @return A new `OncoExperiment` object containing the subset of the original data.
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
    ## ---------------------------------------------------
    ## Validate arguments
    ## ---------------------------------------------------

    if (!is(object, "OncoExperiment")) {
      stop("object must be an OncoExperiment")
    }

    if (!is.logical(subset)) {
      stop("subset must be a logical expression")
    }

    if (!is.logical(select)) {
      stop("select must be a logical expression")
    }

    validObject(object)

    col_data <- MultiAssayExperiment::colData(object)
    primary_ids <- rownames(col_data)

    if (is.null(primary_ids)) {
      stop("Cannot subset because `colData(object)` has no primary IDs", call. = FALSE)
    }

    if (anyDuplicated(primary_ids)) {
      stop("Cannot subset because `colData(object)` has duplicate primary IDs", call. = FALSE)
    }

    if (missing(subset)) {
      # User did not specify a subset, so we keep everything
      keep <- rep(TRUE, length(primary_ids)) # Keep all rows
      warning("No subset expression specified, keeping all rows")
    } else {
      subset_expr <- substitute(subset)
      keep <- .eval_subset(subset_expr, col_data, parent.env()) # Evaluate the subset expression
    }

    keep_ids <- primary_ids[keep]

    if (length(keep_ids) == 0L) {
      warning(
        "Subset expression matched zero models/samples. Returning an empty OncoExperiment.",
        .call = FALSE
      )
    }

    out <- object[, keep_ids, ]
    validObject(out)
    out
  }
)

#' @keywords internal
.eval_subset <- function(
  subset_expr,
  col_data,
  parent_env
) {
  col_data_df <- as.data.frame(col_data, stringsAsFactors = FALSE)

  if (nrow(col_data_df) == 0L) {
    stop("Cannot subset because `colData(object)` is empty", call. = FALSE)
  }

  if (ncol(col_data_df) == 0L) {
    stop("Cannot subset because `colData(object)` has no columns", call. = FALSE)
  }

  # Create a new environment with the column data, prevents namespace conflicts
  subset_env <- list2env(
    as.list(col_data_df),
    parent = parent_env
  )

  keep <- tryCatch(
    eval(subset_expr, envir = subset_env, enclos = parent_env),
    error = function(e) {
      stop("Error in subset expression: ", e$message, call. = FALSE)
    }
  )

  if (!is.logical(keep)) {
    stop("Result of subset expression is not a logical vector", call. = FALSE)
  }

  if (length(keep) != nrow(col_data_df)) {
    stop(
      "Result of subset expression has incorrect length. Expected ",
      nrow(col_data_df), ",
         got ",
      length(keep),
      ".",
      call. = FALSE
    )
  }

  keep
}
