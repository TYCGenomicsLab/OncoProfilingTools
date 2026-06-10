# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

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
    
  }
)