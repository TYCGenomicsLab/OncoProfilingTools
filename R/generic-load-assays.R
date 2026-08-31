# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Load assay data into an OncoExperiment object
#' The function will coerce the assay into the corresponding OncoAssay class
#' before being added to the OncoExperiment object.
#'
#' @param object An OncoExperiment object to which the assay will be added.
#' @param assay_type A character string of the type of assay being loaded. Use the `AssayTypes` enum to easily select a valid assay type.
#' @param data A character string specifying specific file path where the assay data is located.
#' @param metadata A character string specifying the path to the metadata file. This is optional but highly recommended.
#' @return An updated OncoExperiment with loaded assays in the @assays slot.
#' @export
methods::setGeneric(
  name = "load_assays",
  def = function(object,
                 assay_type = NULL,
                 data = NULL,
                 metadata = NULL,
                 ...) {
    if (missing(object)) {
      stop("`object` is required.", call. = FALSE)
    }

    if (missing(assay_type)) {
      stop("`assay_type` is required.", call. = FALSE)
    }

    if (missing(data)) {
      stop("`data` is required.", call. = FALSE)
    }

    if (missing(metadata)) {
      warning("No metadata was provided when calling `load_assays`. Consider providing metadata for a more complete analysis.", call. = FALSE)
    }

    standardGeneric("load_assays")
  }
)
