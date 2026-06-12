# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Load assay data into an OncoExperiment object
#'
#' Loads one or more supported assay files into an `OncoExperiment`.
#'
#' The assay type determines which registered loader is used. Each loader returns
#' an assay object, such as an `ExpressionAssay`, which is then added to the
#' `OncoExperiment`.
#'
#' @param object An `OncoExperiment` object to which the assay will be added.
#' @param assay_type A character vector specifying the assay type(s) being loaded,
#'   such as `"Expression"` or `"PRISM"`.
#' @param path A character string specifying the directory or specific file path
#'   where the assay data is located.
#' @param overwrite A logical value indicating whether to overwrite an existing
#'   assay with the same name.
#' @param ... Additional arguments passed to assay-specific loader functions.
#'
#' @return The updated `OncoExperiment` object.
#'
#' @include registry-assays.R
#' @export
methods::setMethod(
  f = "load_assays",
  signature = c(object = "OncoExperiment"),
  definition = function(
    object,
    assay_type = NULL,
    path = ".cache",
    overwrite = FALSE,
    ...
  ) {
    ## ---------------------------------------------------
    ## Input validation via guard clauses
    ## ---------------------------------------------------

    if (is.null(assay_type)) {
      stop("Please specify an assay type to load.", call. = FALSE)
    }

    if (!is.character(assay_type)) {
      stop("`assay_type` must be a character vector.", call. = FALSE)
    }

    if (length(assay_type) == 0L) {
      stop("`assay_type` cannot be an empty character vector.", call. = FALSE)
    }

    if (anyNA(assay_type) || any(!nzchar(assay_type))) {
      stop("`assay_type` cannot contain NA values or empty strings.", call. = FALSE)
    }

    assay_type <- unique(assay_type)

    if (is.null(path)) {
      stop("Please specify a path to load the assay data from.", call. = FALSE)
    }

    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
      stop("`path` must be a single non-empty character string.", call. = FALSE)
    }

    if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
      stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
    }

    ## ---------------------------------------------------
    ## Validate that the assay type is loadable
    ## ---------------------------------------------------

    loadable_types <- unique(stats::na.omit(.get_loadable_assay_registry()$assay_type))

    unsupported_types <- setdiff(assay_type, loadable_types)

    if (length(unsupported_types) > 0L) {
      stop(
        sprintf(
          "Unsupported or non-loadable assay type(s): %s. Loadable types are: %s",
          paste(unsupported_types, collapse = ", "),
          paste(loadable_types, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    ## ---------------------------------------------------
    ## Resolve files to load based on assay type
    ## ---------------------------------------------------

    if (!file.exists(path)) {
      stop(
        sprintf(
          "Specified path '%s' does not exist. Run `download_assays()` first or provide `path` manually.",
          path
        ),
        call. = FALSE
      )
    }

    resolved_files <- .resolve_files_to_load(
      path = path,
      assay_type = assay_type
    )

    if (!inherits(resolved_files, "data.frame") || nrow(resolved_files) == 0L) {
      stop(
        "No loadable assay files were resolved.",
        call. = FALSE
      )
    }

    ## ---------------------------------------------------
    ## Guard against ambiguous duplicate assay types
    ## ---------------------------------------------------

    duplicated_assay_types <- resolved_files$assay_type[
      duplicated(resolved_files$assay_type)
    ]

    if (length(duplicated_assay_types) > 0L) {
      stop(
        "Multiple loadable files matched the same assay type: ",
        paste(unique(duplicated_assay_types), collapse = ", "),
        "\nMatched files: ",
        paste(resolved_files$filename, collapse = ", "),
        "\nPlease provide a direct file path or use distinct assay types in the registry.",
        call. = FALSE
      )
    }

    ## ---------------------------------------------------
    ## Load each resolved file and add it to object
    ## ---------------------------------------------------

    for (i in seq_len(nrow(resolved_files))) {
      message(sprintf("Loading assay file: %s", resolved_files$filename[[i]]))

      assay <- .load_assay(
        path = resolved_files$path[[i]],
        filename = resolved_files$filename[[i]],
        assay_type = resolved_files$assay_type[[i]],
        OncoAssay_class = resolved_files$OncoAssay_class[[i]],
        loader = resolved_files$loader[[i]],
        ...
      )

      object <- add_assay(
        object = object,
        assay = assay,
        name = resolved_files$assay_type[[i]],
        overwrite = overwrite
      )
    }

    ## ---------------------------------------------------
    ## Final validation and return
    ## ---------------------------------------------------

    validObject(object)

    object
  }
)