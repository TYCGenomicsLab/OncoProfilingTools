# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Add an OncoAssay to an OncoExperiment
#'
#' Adds an object inheriting from `OncoAssay` to the `assays` slot of an
#' `OncoExperiment`.
#'
#' @param object An `OncoExperiment` object.
#' @param assay An object inheriting from `OncoAssay`.
#' @param name Optional character string specifying the assay name. If `NULL`,
#'   `assay@name` is used.
#' @param overwrite Logical value indicating whether to overwrite an existing
#'   assay with the same name.
#' @param ... Additional arguments. Currently unused.
#' @return The updated `OncoExperiment` object.
#' @export
methods::setMethod(
  f = "add_assay",
  signature = c(object = "OncoExperiment", assay = "OncoAssay"),
  definition = function(
    object,
    assay,
    name = NULL,
    overwrite = FALSE,
    ...
  ) {
    ## ---------------------------------------------------
    ## Validate arguments
    ## ---------------------------------------------------

    if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
      stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
    }

    validObject(assay)

    assay_name <- if (is.null(name)) {
      assay@name
    } else {
      name
    }

    if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
      stop(
        "`name` must be NULL or a single non-empty character string.",
        call. = FALSE
      )
    }

    ## ---------------------------------------------------
    ## Validate assay container
    ## ---------------------------------------------------

    if (is.null(object@assays)) {
      object@assays <- list()
    }

    if (!is.list(object@assays)) {
      stop("`object@assays` must be a list.", call. = FALSE)
    }

    ## ---------------------------------------------------
    ## Prevent accidental overwrite
    ## ---------------------------------------------------

    if (assay_name %in% names(object@assays) && !isTRUE(overwrite)) {
      stop(
        "An assay named `",
        assay_name,
        "` already exists in this OncoExperiment. ",
        "Use `overwrite = TRUE` to replace it.",
        call. = FALSE
      )
    }

    ## ---------------------------------------------------
    ## Insert assay
    ## ---------------------------------------------------

    assay@name <- assay_name
    validObject(assay)

    object@assays[[assay_name]] <- assay

    validObject(object)

    object
  }
)
