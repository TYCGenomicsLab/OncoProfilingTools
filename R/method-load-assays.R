# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Load assay data into an OncoExperiment object
#'
#' @description This function loads and parses the specified data and metadata files and appends to the @assays 
#' slot of the `OncoExperiment` object. You can view available assay types by printing `AssayTypes`.
#'
#' @param object An `OncoExperiment` object to which the assay will be added.
#' @param assay_type A character string of the type of assay being loaded. Use the `AssayTypes` enum to easily select a valid assay type.
#' @param data A character string specifying the path to the data file. This is required.
#' @param metadata A character string specifying the path to the metadata file. This is optional but highly recommended.
#' @param overwrite A logical value indicating whether to overwrite an existing
#'   assay with the same name.
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
    data = NULL,
    metadata = NULL,
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

    if (is.null(data)) {
      stop("Please specify a path to load the assay data from.", call. = FALSE)
    }

    if (!is.character(data) || length(data) != 1L || is.na(data) || !nzchar(data)) {
      stop("`data` must be a single non-empty character string.", call. = FALSE)
    }

    if (!is.null(metadata)) {
      if (!is.character(metadata) || length(metadata) != 1L || is.na(metadata) || !nzchar(metadata)) {
        stop("`metadata` must be a single non-empty character string.", call. = FALSE)
      }
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

    if (!file.exists(data)) {
      stop(
        sprintf(
          "Specified data file '%s' does not exist. Please check the file path and try again.",
          data
        ),
        call. = FALSE
      )
    }

    resolved_files <- .resolve_files_to_load(
      path = data,
      assay_type = assay_type
    )

    if (!inherits(resolved_files, "data.frame") || nrow(resolved_files) == 0L) {
      stop(
        "No loadable assay files were resolved.",
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
