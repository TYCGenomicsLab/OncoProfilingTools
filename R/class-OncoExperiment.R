# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' OncoExperiment S4 class definition
#'
#' `OncoExperiment` is a `MultiAssayExperiment` subclass for oncology profiling
#' projects.
#'
#' It stores multiple assay datasets, such as bulk expression, PRISM drug
#' response, mutation data, copy number data, and single-cell experiments, using
#' the standard Bioconductor `MultiAssayExperiment` infrastructure.
#'
#' Additional oncology-specific results, such as fitted models and graphs, are
#' stored in dedicated slots.
#'
#' @slot project_name A character string representing the project name.
#' @slot version A character string representing the package/object version.
#' @slot models A list of fitted models, predictions, or other analytical results.
#' @slot graphs A list of plots or graph objects related to the experiment.
#'
#' @return An object of class `OncoExperiment`.
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
methods::setClass(
  "OncoExperiment",
  contains = "MultiAssayExperiment",
  slots = c(
    project_name = "character",
    version = "character",
    models = "list",
    graphs = "list"
  )
)

#' Create an empty OncoExperiment object
#'
#' Creates an empty oncology-focused `MultiAssayExperiment` subclass.
#'
#' Assays and metadata should be added after construction using package helper
#' functions such as `load_assays()` and `add_assay()`.
#'
#' @param project_name A character string representing the project name.
#'
#' @return An empty object of class `OncoExperiment`.
#'
#' @examples
#' experiment <- OncoExperiment()
#'
#' @name OncoExperiment
#'
#' @export
OncoExperiment <- function(
  project_name = "OncoExperiment Project"
) {
  if (!is.character(project_name) || length(project_name) != 1L || is.na(project_name) || !nzchar(project_name)) {
    stop("`project_name` must be a single non-empty character string.", call. = FALSE)
  }

  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(),
    metadata = list(
      project_name = project_name,
      version = "0.1.0"
    )
  )

  object <- methods::as(mae, "OncoExperiment")

  object@project_name <- project_name
  object@version <- "0.1.0"
  object@models <- list()
  object@graphs <- list()

  validObject(object)

  object
}


#' Show an OncoExperiment object
#'
#' @param object An `OncoExperiment` object.
#' @return Invisibly returns `NULL`.
#' @export
methods::setMethod(
  f = "show",
  signature = signature(object = "OncoExperiment"),
  definition = function(object) {
    cat("An OncoExperiment object\n")
    cat("Project Name:", object@project_name, "\n")
    cat("Version:", object@version, "\n")

    experiments <- MultiAssayExperiment::experiments(object)
    experiment_names <- names(experiments)

    if (is.null(experiment_names) || length(experiment_names) == 0L) {
      cat("Experiments: none\n")
      return(invisible(NULL))
    }

    cat("Experiments:", paste(experiment_names, collapse = ", "), "\n")

    experiment_summary <- lapply(experiment_names, function(experiment_name) {
      experiment <- experiments[[experiment_name]]

      data_dim <- tryCatch(
        dim(experiment),
        error = function(e) NULL
      )

      dim_text <- if (is.null(data_dim)) {
        "unknown"
      } else {
        paste(data_dim, collapse = " x ")
      }

      data.frame(
        experiment = experiment_name,
        class = class(experiment)[1L],
        dim = dim_text,
        stringsAsFactors = FALSE
      )
    })

    experiment_summary <- do.call(rbind, experiment_summary)

    cat("\nExperiment summary:\n")
    print(experiment_summary, row.names = FALSE)

    if (length(object@models) > 0L) {
      cat("\nModels:", paste(names(object@models), collapse = ", "), "\n")
    }

    if (length(object@graphs) > 0L) {
      cat("Graphs:", paste(names(object@graphs), collapse = ", "), "\n")
    }

    invisible(NULL)
  }
)
