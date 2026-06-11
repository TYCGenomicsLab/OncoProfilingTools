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
methods::setClass("ProteinExpressionAssay",
  contains = "OncoAssay",
  slots = c(
    unit = "character",
    normalized = "logical"
  )
)

#' @keywords internal
#' Validity method for narrowing data slot to a numeric matrix with row and column names.
methods::setValidity("ProteinExpressionAssay", function(object) {
  if (!is.matrix(object@data)) {
    return("ExpressionAssay `data` must be a matrix.")
  }

  if (!is.numeric(object@data)) {
    return("ExpressionAssay `data` must be a numeric matrix.")
  }

  if (is.null(rownames(object@data))) {
    return("ExpressionAssay `data` must have row names for samples/models.")
  }

  if (is.null(colnames(object@data))) {
    return("ExpressionAssay `data` must have column names for genes/features.")
  }

  TRUE
})

#' @export
ProteinExpressionAssay <- function(
  data, metadata = list(), row_data = list(),
  col_data = list(), layers = list(), name = "ExpressionAssay",
  unit = "TPM", normalized = TRUE
) {
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  new("ProteinExpressionAssay",
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

#' Show an ExpressionAssay object
#'
#' @param object An `ExpressionAssay` object.
#'
#' @return Invisibly returns `NULL`.
#' @export
methods::setMethod(
  f = "show",
  signature = signature(object = "ProteinExpressionAssay"),
  definition = function(object) {
    cat("An ExpressionAssay object\n")
    cat("Name:", object@name, "\n")
    cat("Unit:", object@unit, "\n")
    cat("Normalized:", object@normalized, "\n")

    data_dim <- dim(object@data)

    cat("Data:", paste(data_dim, collapse = " x "), class(object@data)[1L], "\n")

    if (!is.null(rownames(object@data))) {
      cat("Rows:", length(rownames(object@data)), "models/samples\n")
    }

    if (!is.null(colnames(object@data))) {
      cat("Columns:", length(colnames(object@data)), "features\n")
    }

    if (length(object@metadata) > 0L) {
      cat("Metadata fields:", paste(names(object@metadata), collapse = ", "), "\n")
    }

    invisible(NULL)
  }
)
