# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' List assay names
#'
#' Returns the names of assays stored in an object.
#'
#' @param object An object containing assays.
#' @param ... Additional arguments passed to methods.
#'
#' @return A character vector of assay names.
#' @include class-OncoExperiment.R
#' @export
methods::setGeneric(
  name = "assays",
  def = function(object, ...) {
    standardGeneric("assays")
  }
)

#' Extract an assay
#'
#' Extracts a single assay from an object by name.
#'
#' @param object An object containing assays.
#' @param name A single character string specifying the assay name.
#' @param ... Additional arguments passed to methods.
#'
#' @return An assay object.
#' @include class-OncoExperiment.R
#' @export
methods::setGeneric(
  name = "assay",
  def = function(object, name, ...) {
    standardGeneric("assay")
  }
)

#' List assay names in an OncoExperiment
#'
#' Returns the names of assays stored in the `assays` slot of an
#' `OncoExperiment`.
#'
#' @param object An `OncoExperiment` object.
#' @param ... Additional arguments. Currently unused.
#'
#' @return A character vector of assay names.
#'
#' @examples
#' \dontrun{
#' assays(object)
#' }
#'
#' @export
methods::setMethod(
  f = "assays",
  signature = signature(object = "OncoExperiment"),
  definition = function(object, ...) {
    names(object@assays)
  }
)

#' Extract an assay from an OncoExperiment
#'
#' Extracts a single assay from the `assays` slot of an `OncoExperiment` by
#' name.
#'
#' @param object An `OncoExperiment` object.
#' @param name A single character string specifying the name of the assay to
#'   extract.
#' @param ... Additional arguments. Currently unused.
#'
#' @return An object inheriting from `OncoAssay`.
#'
#' @examples
#' \dontrun{
#' expr <- assay(object, "Expression")
#' }
#'
#' @export
methods::setMethod(
  f = "assay",
  signature = signature(object = "OncoExperiment"),
  definition = function(object, name, ...) {
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
      stop("`name` must be a single non-empty character string.", call. = FALSE)
    }

    if (!name %in% names(object@assays)) {
      stop(
        "Assay `",
        name,
        "` not found. Available assays: ",
        paste(names(object@assays), collapse = ", "),
        call. = FALSE
      )
    }

    object@assays[[name]]
  }
)
