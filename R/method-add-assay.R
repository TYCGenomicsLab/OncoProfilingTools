#' Add an assay experiment to an OncoExperiment
#'
#' Adds a named assay experiment to the `experiments()` of an `OncoExperiment`.
#'
#' @param object An `OncoExperiment` object.
#' @param assay An assay object, such as a `SummarizedExperiment`,
#'   `ExpressionAssay`, or `SingleCellExperiment`.
#' @param name A single character string specifying the experiment name.
#' @param overwrite Logical value indicating whether to overwrite an existing
#'   experiment with the same name.
#' @param ... Additional arguments. Currently unused.
#'
#' @return The updated `OncoExperiment` object.
#' @export
methods::setMethod(
  f = "add_assay",
  signature = c(object = "OncoExperiment", assay = "SummarizedExperiment"),
  definition = function(
    object,
    assay,
    name = NULL,
    overwrite = FALSE,
    ...
  ) {
    ## ---------------------------------------------------
    ## Validate arguments
    ## ---------------------------------------------------

    if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
      stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
    }

    if (is.null(name) || !is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
      stop("`name` must be a single non-empty character string.", call. = FALSE)
    }

    validObject(assay)

    ## ---------------------------------------------------
    ## Get existing experiments
    ## ---------------------------------------------------

    experiments <- MultiAssayExperiment::experiments(object)

    if (!inherits(experiments, "ExperimentList")) {
      experiments <- MultiAssayExperiment::ExperimentList(experiments)
    }

    if (name %in% names(experiments) && !isTRUE(overwrite)) {
      stop(
        "An assay named `",
        name,
        "` already exists in this OncoExperiment. ",
        "Use `overwrite = TRUE` to replace it.",
        call. = FALSE
      )
    }

    experiments[[name]] <- assay

    ## ---------------------------------------------------
    ## Rebuild MultiAssayExperiment
    ## ---------------------------------------------------
    ## This lets MultiAssayExperiment regenerate sampleMap correctly from the
    ## experiment colnames instead of trying to patch the old empty sampleMap.

    mae <- MultiAssayExperiment::MultiAssayExperiment(
      experiments = experiments,
      metadata = S4Vectors::metadata(object)
    )

    updated <- methods::as(mae, "OncoExperiment")

    ## ---------------------------------------------------
    ## Preserve OncoExperiment-specific slots
    ## ---------------------------------------------------

    updated@project_name <- object@project_name
    updated@version <- object@version
    updated@models <- object@models
    updated@graphs <- object@graphs

    validObject(updated)

    updated
  }
)