# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Add an assay to an object
#'
#' Adds an assay object to a container object.
#'
#' @param object The object to which the assay will be added.
#' @param assay The assay object to add.
#' @param name Optional character string specifying the assay name. If `NULL`,
#'   the assay object's existing `name` slot is used.
#' @param overwrite Logical value indicating whether to overwrite an existing
#'   assay with the same name.
#' @param ... Additional arguments passed to methods.
#'
#' @return The updated object.
#'
#' @export
methods::setGeneric(
  name = "add_assay",
  def = function(object, assay, name = NULL, overwrite = FALSE, ...) {
    standardGeneric("add_assay")
  }
)