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

    # Keep any top-level metadata already attached to the OncoExperiment so
    # rebuilding the container does not drop user-facing columns.
    existing_col_data <- MultiAssayExperiment::colData(object)

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

    # Collect metadata from the container and each assay, then merge them into
    # a single aligned colData table for the rebuilt MultiAssayExperiment.
    col_data_sources <- c(
      list(existing_col_data),
      lapply(experiments, function(experiment) {
        MultiAssayExperiment::colData(experiment)
      })
    )

    primary_ids <- unique(unlist(lapply(col_data_sources, rownames), use.names = FALSE))
    primary_ids <- primary_ids[!is.na(primary_ids) & nzchar(primary_ids)]

    merged_col_data <- S4Vectors::DataFrame(row.names = primary_ids)

    for (source_col_data in col_data_sources) {
      source_col_data <- S4Vectors::DataFrame(source_col_data)

      if (nrow(source_col_data) == 0L || ncol(source_col_data) == 0L) {
        next
      }

      source_row_names <- rownames(source_col_data)

      if (is.null(source_row_names)) {
        next
      }

      for (col_name in colnames(source_col_data)) {
        # Align each metadata column to the final primary IDs so the rows stay
        # matched to the correct models after reconstruction.
        aligned_values <- source_col_data[[col_name]][match(primary_ids, source_row_names)]

        if (!col_name %in% colnames(merged_col_data)) {
          merged_col_data[[col_name]] <- aligned_values
        } else {
          existing_values <- merged_col_data[[col_name]]
          replace_idx <- is.na(existing_values) & !is.na(aligned_values)
          existing_values[replace_idx] <- aligned_values[replace_idx]
          merged_col_data[[col_name]] <- existing_values
        }
      }
    }

    mae <- MultiAssayExperiment::MultiAssayExperiment(
      experiments = experiments,
      colData = merged_col_data,
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
