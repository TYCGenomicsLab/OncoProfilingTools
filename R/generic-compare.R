# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

methods::setGeneric(
  name = "compare",
  def = function(object, group1, group2, unit, p_adj_method = "BH", effect_threshold = 0, ...),
  signature = c(
    object = "OncoExperiment",
    group1 = "OncoExperiment",
    group2 = "OncoExperiment",
    unit = "character",
    p_adj_method = "character",
    effect_threshold = "numeric"
  ) {
    standardGeneric("compare")
  }
)