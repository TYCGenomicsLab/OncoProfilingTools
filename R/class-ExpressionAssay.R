# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' ExpressionAssay S4 class definition
#' This class is designed to hold gene expression assay data, including the expression matrix, metadata, and additional information relevant to the assay.
#' It inherits from the OncoAssay virtual class. Usually, you won't create an ExpressionAssay manually but use a helper function for your data source.
#'
#' Typical use case is for DepMap gene expression data.
#' @include class-OncoAssay.R
#' @name ExpressionAssay
#' @slot unit A character string representing the unit of measurement for the expression data, such as "TPM", "RPKM", or "raw counts".
#' @slot normalized A logical value indicating whether the expression data has been normalized (TRUE) or is in raw form (FALSE).
#' @return An object of class ExpressionAssay
#' @inherit OncoAssay
#' @export
methods::setClass("ExpressionAssay",
  contains = "OncoAssay",
  slots = c(
    unit = "character",
    normalized = "logical"
  )
)

#' @export
ExpressionAssay <- function(data, metadata = list(), row_data = list(), col_data = list(), layers = list(), name = "ExpressionAssay", unit = "TPM", normalized = TRUE) {
  new("ExpressionAssay",
    data = data,
    metadata = metadata,
    row_data = row_data,
    col_data = col_data,
    layers = layers,
    name = name,
    unit = unit,
    normalized = normalized
  )
}
