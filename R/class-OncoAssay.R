# Author: Jason LaPierre
# Last update: June 4, 2026
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' OncoAssay virtual S4 class definition
#' This class serves as a base class for different types of assays in oncology experiments,
#' such as gene expression matrices, mutation data, or other experimental data. It is designed
#' to hold the assay data, metadata, and any additional information relevant to the assay.
#' @name OncoAssay
#' @slot metadata A list of metadata information related to the assay, such as sample annotations, experimental conditions, or other relevant information.
#' @slot data A matrix containing the assay data, such as gene expression values or mutation counts.
#' @slot row_data A list of metadata information for the rows of the data matrix, such as gene annotations or mutation information.
#' @slot col_data A list of metadata information for the columns of the data matrix, such as sample annotations or experimental conditions.
#' @slot layers A list of additional data layers related to the assay, such as normalized data, log-transformed data, or other derived data matrices.
#' @slot name A character string representing the name of the assay.
#' @return An object of class OncoAssay
#' @export
methods::setClass("OncoAssay",
  contains = "VIRTUAL",
  slots = c(
    metadata = "list",
    data = "matrix",
    row_data = "list",
    col_data = "list",
    layers = "list",
    name = "character"
  )
)
