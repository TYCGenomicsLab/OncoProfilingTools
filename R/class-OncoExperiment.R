# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' OncoExperiment S4 class definition
#' This class is designed to hold all relevant data and results for an oncology experiment,
#' including assays, models, graphs, and metadata.
#' @slot Optional project_name A character string representing the name of the project.
#' @slot version A character string representing the version of the package.
#' @slot assays A list of OncoAssay objects, such as gene expression matrices or other experimental data.
#' @slot models A list of OncoModel objects, which may include fitted models, predictions, or other analytical results.
#' @slot graphs A list of OncoGraph objects, which may include ggplot2 objects or other visualizations related to the experiment.
#' @slot metadata A list of metadata information, such as sample annotations, experimental conditions, or other relevant information.
#' @return An object of class OncoExperiment
#' @examples
#' # Create an empty OncoExperiment object
#' experiment <- OncoExperiment(project_name = "My Oncology Experiment")
#' @export
methods::setClass("OncoExperiment",
  slots = c(
    project_name = "character",
    version = "character",
    assays = "list",
    models = "list",
    graphs = "list",
    metadata = "list"
  )
)

#' OncoExperiment constructor function
#' This function creates a new OncoExperiment object with the specified project name, assays,
#' models, graphs, and metadata. The version is set to "0.1.0" by default.
#' This class is designed to hold all relevant data and results for an oncology experiment,
#' including assays, models, graphs, and metadata. You should not need to define any slots for a blank experiment.
#' @param project_name A character string representing the name of the project. Default is "OncoExperiment Project".
#' @param assays A list of OncoAssay objects, such as gene expression matrices or other experimental data. Default is an empty list.
#' @param models A list of OncoModel objects, which may include fitted models, predictions, or other analytical results. Default is an empty list.
#' @param graphs A list of OncoGraph objects, which may include ggplot2 objects or other visualizations related to the experiment. Default is an empty list.
#' @param metadata A list of metadata information, such as sample annotations, experimental conditions, or other relevant information. Default is an empty list.
#' @return An object of class OncoExperiment
#' @examples
#' # Create an empty OncoExperiment object with default values
#' experiment <- OncoExperiment()
#' @export
OncoExperiment <- function(
  project_name = "OncoExperiment Project",
  assays = list(),
  models = list(),
  graphs = list(),
  metadata = list()
) {
  # TODO: add validation checks for input parameters once sub-classes are defined

  methods::new("OncoExperiment",
    project_name = project_name,
    version = "0.1.0",
    assays = assays,
    models = models,
    graphs = graphs,
    metadata = metadata
  )
}

#' Show method for OncoExperiment objects.
#' This method defines how OncoExperiment objects are displayed when printed.
setMethod("show", "OncoExperiment", function(object) {
  cat("An OncoExperiment object\n")
  cat("Project Name:", object@project_name, "\n")
  cat("Version:", object@version, "\n")
})
