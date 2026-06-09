# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Load assay data into an OncoExperiment object
#' The function will coerce the assay into the corresponding OncoAssay class
#' before being added to the OncoExperiment object.
#'
#' @param object An OncoExperiment object to which the assay will be added.
#' @param assay_type A character string specifying the type of assay being loaded (e.g., "Expression", "PRISM", etc.). This will determine the specific OncoAssay used and correct loading pattern.
#' @param path A character string specifying the directory or specific file path where the assay data is located. The function will look for files in this location that match the expected naming conventions for the specified assay type.
#' @return An updated OncoExperiment with loaded assays in the @assays slot.
#' @export
methods::setGeneric(
  name="load_assays",
  def=function(object, assay_type=NULL, path=NULL, ...) {
    if (missing(object)) {
      stop("`object` is required.", call. = FALSE)
    }

    standardGeneric("load_assays")
  }
)