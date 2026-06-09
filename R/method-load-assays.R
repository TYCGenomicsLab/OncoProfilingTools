# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Load assay data into an OncoExperiment object
#'
#' The function will coerce the assay into the corresponding OncoAssay class
#' before being added to the OncoExperiment object. The assay type specified will determine the specific OncoAssay used and the correct loading pattern.
#'
#' @param object An OncoExperiment object to which the assay will be added.
#' @param assay_type A character string specifying the type of assay being loaded (e.g., "Expression", "PRISM", etc.). This will determine the specific OncoAssay used and correct loading pattern.
#' @param path A character string specifying the directory or specific file path where the assay data is located. The function will look for files in this location that match the expected naming conventions for the specified assay type.
#' @param overwrite A logical value indicating whether to overwrite an existing assay of the same type in the OncoExperiment object. Default is FALSE, which will prevent overwriting and throw an error if an assay of the same type already exists.
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
      stop("Assay type must be a character string.", call. = FALSE)
    }

    if (length(assay_type) == 0L) {
      stop("Assay type cannot be an empty string.", call. = FALSE)
    }

    if (is.null(path)) {
      stop("Please specify a path to load the assay data from.", call. = FALSE)
    }

    ## ---------------------------------------------------
    ## Validate that the assay type is supported
    ## ---------------------------------------------------
    supported_types <- get_supported_types()
    if (!(assay_type %in% supported_types)) {
      stop(
        sprintf(
          "Unsupported assay type '%s'. Supported types are: %s",
          assay_type,
          paste(supported_types, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    ## ---------------------------------------------------
    ## Resolve files to load based on assay type
    ## ---------------------------------------------------
    if (!dir.exists(path) | !file.exists(path)) {
      stop(
        sprintf(
          "Specified path '%s' does not exist. Run `download_assays()` first or provide `path` manually.",
          search_path
        ),
        call. = FALSE
      )
    }


    files_to_load <- .resolve_files_to_load(path, assay_type)

    message(files_to_load)
  }
)
