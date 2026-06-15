# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

#' Find a metadata sidecar file beside an assay file
#'
#' @param path Path to the primary assay file.
#' @param filename Metadata filename to find.
#'
#' @return A normalized character path if found, otherwise `NULL`.
#'
#' @keywords internal
.find_sidecar_file <- function(path, filename) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(filename) || length(filename) != 1L || is.na(filename) || !nzchar(filename)) {
    stop("`filename` must be a single non-empty character string.", call. = FALSE)
  }

  candidate <- file.path(dirname(path), filename)

  if (file.exists(candidate)) {
    normalizePath(candidate, winslash = "/", mustWork = TRUE)
  } else {
    NULL
  }
}

#' Load DepMap Model.csv metadata
#'
#' Loads and aligns DepMap `Model.csv` metadata to assay column names.
#'
#' @param model_metadata_path Path to `Model.csv`, or `NULL`.
#' @param model_ids Character vector of model IDs, usually assay column names.
#'
#' @return An `S4Vectors::DataFrame` aligned to `model_ids`.
#'
#' @keywords internal
.load_model_metadata <- function(
  model_metadata_path = NULL,
  model_ids
) {
  if (!is.character(model_ids) || length(model_ids) == 0L) {
    stop("`model_ids` must be a non-empty character vector.", call. = FALSE)
  }

  if (anyNA(model_ids) || any(!nzchar(model_ids))) {
    stop("`model_ids` cannot contain NA values or empty strings.", call. = FALSE)
  }

  if (anyDuplicated(model_ids)) {
    duplicated_models <- unique(model_ids[duplicated(model_ids)])

    stop(
      "`model_ids` contains duplicated values. Examples: ",
      paste(utils::head(duplicated_models, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(model_metadata_path)) {
    return(S4Vectors::DataFrame(
      ModelID = model_ids,
      row.names = model_ids
    ))
  }

  if (!is.character(model_metadata_path) || length(model_metadata_path) != 1L || is.na(model_metadata_path) || !nzchar(model_metadata_path)) {
    stop("`model_metadata_path` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(model_metadata_path)) {
    stop("Model metadata file does not exist: ", model_metadata_path, call. = FALSE)
  }

  message("Reading model metadata file: ", basename(model_metadata_path))

  model_df <- data.table::fread(
    file = model_metadata_path,
    check.names = FALSE,
    data.table = FALSE,
    showProgress = interactive()
  )

  if (nrow(model_df) == 0L) {
    stop("Model metadata file is empty: ", model_metadata_path, call. = FALSE)
  }

  if (!"ModelID" %in% colnames(model_df)) {
    stop("Model metadata file must contain a `ModelID` column.", call. = FALSE)
  }

  if (anyNA(model_df$ModelID) || any(!nzchar(model_df$ModelID))) {
    stop("`Model.csv` contains NA or empty `ModelID` values.", call. = FALSE)
  }

  if (anyDuplicated(model_df$ModelID)) {
    duplicated_models <- unique(model_df$ModelID[duplicated(model_df$ModelID)])

    stop(
      "`Model.csv` contains duplicated ModelID values. Examples: ",
      paste(utils::head(duplicated_models, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  match_idx <- match(model_ids, model_df$ModelID)

  if (anyNA(match_idx)) {
    missing_models <- model_ids[is.na(match_idx)]

    warning(
      "Some assay model IDs were not found in Model.csv. Missing examples: ",
      paste(utils::head(missing_models, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  model_metadata <- model_df[match_idx, , drop = FALSE]

  ## Preserve exact assay column order and all model IDs.
  model_metadata$ModelID <- model_ids
  rownames(model_metadata) <- model_ids

  S4Vectors::DataFrame(model_metadata)
}
