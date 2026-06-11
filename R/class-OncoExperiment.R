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

#' Show an OncoExperiment object
#'
#' @param object An `OncoExperiment` object.
#'
#' @return Invisibly returns `NULL`.
#' @export
methods::setMethod(
  f = "show",
  signature = signature(object = "OncoExperiment"),
  definition = function(object) {
    cat("An OncoExperiment object\n")

    if ("project_name" %in% slotNames(object)) {
      cat("Project Name:", object@project_name, "\n")
    }

    if ("version" %in% slotNames(object)) {
      cat("Version:", object@version, "\n")
    }

    assay_names <- names(object@assays)

    if (is.null(assay_names) || length(assay_names) == 0L) {
      cat("Assays: none\n")
      return(invisible(NULL))
    }

    cat("Assays:", paste(assay_names, collapse = ", "), "\n")

    assay_summary <- lapply(assay_names, function(assay_name) {
      assay <- object@assays[[assay_name]]

      data_dim <- tryCatch(
        dim(assay@data),
        error = function(e) NULL
      )

      dim_text <- if (is.null(data_dim)) {
        "unknown"
      } else {
        paste(data_dim, collapse = " x ")
      }

      data_class <- class(assay@data)[1L]

      data.frame(
        assay = assay_name,
        class = class(assay)[1L],
        data_class = data_class,
        dim = dim_text,
        stringsAsFactors = FALSE
      )
    })

    assay_summary <- do.call(rbind, assay_summary)

    cat("\nAssay summary:\n")
    print(assay_summary, row.names = FALSE)

    invisible(NULL)
  }
)
