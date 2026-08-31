# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab


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

#' Load PRISM cell line metadata
#'
#' Loads and aligns PRISM cell line metadata to drug response assay column names.
#'
#' The PRISM cell line metadata is expected to contain a `depmap_id` column
#' matching the model IDs in the drug response matrix columns. Because PRISM
#' metadata may contain multiple rows per `depmap_id` across pools, cultures,
#' and screens, duplicated `depmap_id` rows are collapsed to one row per model.
#'
#' @param cell_line_metadata_path Path to the PRISM cell line metadata file, or
#'   `NULL`.
#' @param model_ids Character vector of model IDs, usually assay column names.
#'
#' @return An `S4Vectors::DataFrame` aligned to `model_ids`.
#'
#' @keywords internal
.load_prism_cell_line_metadata <- function(
  cell_line_metadata_path = NULL,
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

  if (is.null(cell_line_metadata_path)) {
    return(S4Vectors::DataFrame(
      model_id = model_ids,
      row.names = model_ids
    ))
  }

  if (!is.character(cell_line_metadata_path) || length(cell_line_metadata_path) != 1L || is.na(cell_line_metadata_path) || !nzchar(cell_line_metadata_path)) {
    stop(
      "`cell_line_metadata_path` must be NULL or a single non-empty character string.",
      call. = FALSE
    )
  }

  if (!file.exists(cell_line_metadata_path)) {
    stop(
      "PRISM cell line metadata file does not exist: ",
      cell_line_metadata_path,
      call. = FALSE
    )
  }

  message("Reading PRISM cell line metadata file: ", basename(cell_line_metadata_path))

  cell_line_df <- data.table::fread(
    file = cell_line_metadata_path,
    check.names = FALSE,
    data.table = FALSE,
    showProgress = interactive()
  )

  if (nrow(cell_line_df) == 0L) {
    stop(
      "PRISM cell line metadata file is empty: ",
      cell_line_metadata_path,
      call. = FALSE
    )
  }

  if (!"depmap_id" %in% colnames(cell_line_df)) {
    stop(
      "PRISM cell line metadata must contain a `depmap_id` column.",
      call. = FALSE
    )
  }

  if (anyNA(cell_line_df$depmap_id) || any(!nzchar(cell_line_df$depmap_id))) {
    stop(
      "PRISM cell line metadata contains NA or empty `depmap_id` values.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Collapse duplicated depmap_id rows
  ## ---------------------------------------------------

  n_rows_before_dedup <- nrow(cell_line_df)

  cell_line_df <- unique(cell_line_df)

  n_exact_duplicates_removed <- n_rows_before_dedup - nrow(cell_line_df)

  collapse_values <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    x <- unique(x)

    if (length(x) == 0L) {
      NA_character_
    } else {
      paste(x, collapse = "; ")
    }
  }

  if (anyDuplicated(cell_line_df$depmap_id)) {
    duplicated_model_ids <- unique(cell_line_df$depmap_id[duplicated(cell_line_df$depmap_id)])

    warning(
      "PRISM cell line metadata contains duplicated `depmap_id` values. ",
      "Collapsing duplicated metadata rows to one row per `depmap_id`. ",
      "Example duplicated IDs: ",
      paste(utils::head(duplicated_model_ids, 10L), collapse = ", "),
      call. = FALSE
    )

    metadata_cols <- setdiff(colnames(cell_line_df), "depmap_id")

    cell_line_df <- cell_line_df |>
      dplyr::group_by(.data$depmap_id) |>
      dplyr::summarise(
        n_metadata_rows_collapsed = dplyr::n(),
        dplyr::across(
          dplyr::all_of(metadata_cols),
          collapse_values
        ),
        .groups = "drop"
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
  } else {
    cell_line_df$n_metadata_rows_collapsed <- 1L
  }

  if (anyDuplicated(cell_line_df$depmap_id)) {
    stop(
      "Internal error: PRISM cell line metadata still contains duplicated `depmap_id` values after collapsing.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Align metadata to assay matrix columns
  ## ---------------------------------------------------

  match_idx <- match(model_ids, cell_line_df$depmap_id)

  if (anyNA(match_idx)) {
    missing_models <- model_ids[is.na(match_idx)]

    warning(
      "Some PRISM model IDs were not found in cell line metadata. Missing examples: ",
      paste(utils::head(missing_models, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  model_metadata <- cell_line_df[match_idx, , drop = FALSE]

  ## Preserve exact assay column order.
  model_metadata$model_id <- model_ids
  rownames(model_metadata) <- model_ids

  ## Add stable aliases while preserving original PRISM columns.
  if ("depmap_id" %in% colnames(model_metadata)) {
    model_metadata$depmap_id <- model_metadata$depmap_id
  }

  if ("ccle_name" %in% colnames(model_metadata)) {
    model_metadata$ccle_name <- model_metadata$ccle_name
  }

  if ("row_id" %in% colnames(model_metadata)) {
    model_metadata$prism_row_id <- model_metadata$row_id
  }

  key_cols <- c(
    "model_id",
    "depmap_id",
    "ccle_name",
    "prism_row_id",
    "pool_id",
    "culture",
    "screen",
    "n_metadata_rows_collapsed"
  )

  other_cols <- setdiff(colnames(model_metadata), key_cols)

  model_metadata <- model_metadata[
    ,
    c(intersect(key_cols, colnames(model_metadata)), other_cols),
    drop = FALSE
  ]

  model_metadata <- S4Vectors::DataFrame(model_metadata)

  S4Vectors::metadata(model_metadata)$n_exact_duplicates_removed <- n_exact_duplicates_removed
  S4Vectors::metadata(model_metadata)$cell_line_metadata_path <- normalizePath(
    cell_line_metadata_path,
    winslash = "/",
    mustWork = TRUE
  )

  model_metadata
}

#' Load PRISM compound metadata
#'
#' Loads and aligns PRISM compound metadata to drug response assay row names.
#'
#' The PRISM compound list is expected to contain an `IDs` column matching the
#' treatment IDs in the response matrix. If the compound metadata contains
#' duplicated `IDs`, rows are collapsed to one row per `IDs` by concatenating
#' unique non-missing values with `"; "`.
#'
#' @param compound_metadata_path Path to the PRISM compound metadata file, or
#'   `NULL`.
#' @param treatment_ids Character vector of treatment IDs, usually assay row
#'   names.
#'
#' @return An `S4Vectors::DataFrame` aligned to `treatment_ids`.
#'
#' @keywords internal
.load_prism_compound_metadata <- function(
  compound_metadata_path = NULL,
  treatment_ids
) {
  if (!is.character(treatment_ids) || length(treatment_ids) == 0L) {
    stop("`treatment_ids` must be a non-empty character vector.", call. = FALSE)
  }

  if (anyNA(treatment_ids) || any(!nzchar(treatment_ids))) {
    stop("`treatment_ids` cannot contain NA values or empty strings.", call. = FALSE)
  }

  if (anyDuplicated(treatment_ids)) {
    duplicated_treatments <- unique(treatment_ids[duplicated(treatment_ids)])

    stop(
      "`treatment_ids` contains duplicated values. Examples: ",
      paste(utils::head(duplicated_treatments, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(compound_metadata_path)) {
    return(S4Vectors::DataFrame(
      treatment_id = treatment_ids,
      row.names = treatment_ids
    ))
  }

  if (!is.character(compound_metadata_path) || length(compound_metadata_path) != 1L || is.na(compound_metadata_path) || !nzchar(compound_metadata_path)) {
    stop(
      "`compound_metadata_path` must be NULL or a single non-empty character string.",
      call. = FALSE
    )
  }

  if (!file.exists(compound_metadata_path)) {
    stop(
      "Compound metadata file does not exist: ",
      compound_metadata_path,
      call. = FALSE
    )
  }

  message("Reading PRISM compound metadata file: ", basename(compound_metadata_path))

  compound_df <- data.table::fread(
    file = compound_metadata_path,
    check.names = FALSE,
    data.table = FALSE,
    showProgress = interactive()
  )

  if (nrow(compound_df) == 0L) {
    stop("Compound metadata file is empty: ", compound_metadata_path, call. = FALSE)
  }

  if (!"IDs" %in% colnames(compound_df)) {
    stop("PRISM compound metadata file must contain an `IDs` column.", call. = FALSE)
  }

  compound_df$IDs <- as.character(compound_df$IDs)

  if (anyNA(compound_df$IDs) || any(!nzchar(compound_df$IDs))) {
    stop("PRISM compound metadata contains NA or empty `IDs` values.", call. = FALSE)
  }

  ## ---------------------------------------------------
  ## Collapse duplicated IDs to one metadata row per ID
  ## ---------------------------------------------------

  n_rows_before_dedup <- nrow(compound_df)

  compound_df <- unique(compound_df)

  n_exact_duplicates_removed <- n_rows_before_dedup - nrow(compound_df)

  collapse_values <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    x <- unique(x)

    if (length(x) == 0L) {
      NA_character_
    } else {
      paste(x, collapse = "; ")
    }
  }

  if (anyDuplicated(compound_df$IDs)) {
    duplicated_ids <- unique(compound_df$IDs[duplicated(compound_df$IDs)])

    warning(
      "PRISM compound metadata contains duplicated `IDs`. ",
      "Collapsing duplicated metadata rows to one row per `IDs`. ",
      "Example duplicated IDs: ",
      paste(utils::head(duplicated_ids, 10L), collapse = ", "),
      call. = FALSE
    )

    metadata_cols <- setdiff(colnames(compound_df), "IDs")

    compound_df <- compound_df |>
      dplyr::group_by(.data$IDs) |>
      dplyr::summarise(
        n_metadata_rows_collapsed = dplyr::n(),
        dplyr::across(
          dplyr::all_of(metadata_cols),
          collapse_values
        ),
        .groups = "drop"
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
  } else {
    compound_df$n_metadata_rows_collapsed <- 1L
  }

  if (anyDuplicated(compound_df$IDs)) {
    stop(
      "Internal error: PRISM compound metadata still contains duplicated `IDs` after collapsing.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Align metadata to response matrix treatment IDs
  ## ---------------------------------------------------

  match_idx <- match(treatment_ids, compound_df$IDs)

  if (anyNA(match_idx)) {
    missing_treatments <- treatment_ids[is.na(match_idx)]

    warning(
      "Some treatment IDs were not found in PRISM compound metadata. Missing examples: ",
      paste(utils::head(missing_treatments, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  compound_metadata <- compound_df[match_idx, , drop = FALSE]

  ## Preserve exact assay row order and all treatment IDs.
  compound_metadata$treatment_id <- treatment_ids
  rownames(compound_metadata) <- treatment_ids

  ## Add stable aliases while preserving original PRISM columns.
  if ("IDs" %in% colnames(compound_metadata)) {
    compound_metadata$compound_id <- compound_metadata$IDs
  }

  if ("Drug.Name" %in% colnames(compound_metadata)) {
    compound_metadata$drug_name <- compound_metadata$Drug.Name
  }

  if ("MOA" %in% colnames(compound_metadata)) {
    compound_metadata$moa <- compound_metadata$MOA
  }

  if ("repurposing_target" %in% colnames(compound_metadata)) {
    compound_metadata$target <- compound_metadata$repurposing_target
  }

  ## Keep important standardized columns first.
  key_cols <- c(
    "treatment_id",
    "compound_id",
    "drug_name",
    "moa",
    "target",
    "IDs",
    "n_metadata_rows_collapsed"
  )

  selected_cols <- c(
    intersect(key_cols, colnames(compound_metadata)),
    setdiff(colnames(compound_metadata), key_cols)
  )

  compound_metadata <- compound_metadata[
    ,
    selected_cols,
    drop = FALSE
  ]

  compound_metadata <- S4Vectors::DataFrame(compound_metadata)

  S4Vectors::metadata(compound_metadata)$n_exact_duplicates_removed <- n_exact_duplicates_removed
  S4Vectors::metadata(compound_metadata)$compound_metadata_path <- normalizePath(
    compound_metadata_path,
    winslash = "/",
    mustWork = TRUE
  )

  compound_metadata
}
