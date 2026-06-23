# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Subset an OncoExperiment
#'
#' @export
methods::setGeneric(
  name = "subset",
  def = function(object,
                 subset,
                 select = TRUE,
                 ...) {
    standardGeneric("subset")
  }
)
