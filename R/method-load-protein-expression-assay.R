# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(data.table)

#' Load a DepMap protein-coding gene expression assay
#'
#' Loads a DepMap protein-coding gene expression file into a
#' `ProteinExpressionAssay`.
#'
#' The DepMap expression file is read, filtered to default model entries, and
#' converted into a Bioconductor-style expression matrix with features as rows
#' and models as columns. If `Model.csv` is present beside the assay file, it is
#' loaded and aligned into `colData()`.
#'
#' @param path A character string specifying the path to the expression assay file.
#' @param assay_name A character string specifying the assay name used in metadata.
#' @param model_metadata_path Optional path to `Model.csv`. If `NULL`, the loader
#'   searches for `Model.csv` beside `path`.
#' @param unit A character string specifying the expression unit.
#' @param normalized A logical value indicating whether the expression data are normalized.
#' @param id_col A character string specifying the expression file column to use
#'   as model/sample identifiers.
#' @param feature_type A character string describing the expression feature type.
#'
#' @return A `ProteinExpressionAssay`, which inherits from
#'   `SummarizedExperiment`.
#'
#' @include internal-metadata.R
#' @keywords internal
load_protein_expression_assay <- function(
  path,
  assay_name = "Protein Expression",
  model_metadata_path = NULL,
  unit = "log2(TPM+1)",
  normalized = TRUE,
  id_col = "ModelID",
  feature_type = "protein_coding_gene"
) {
  ## ---------------------------------------------------
  ## Validate arguments
  ## ---------------------------------------------------

  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop("Protein expression assay file does not exist: ", path, call. = FALSE)
  }

  if (!is.null(model_metadata_path)) {
    if (!is.character(model_metadata_path) || length(model_metadata_path) != 1L || is.na(model_metadata_path) || !nzchar(model_metadata_path)) {
      stop("`model_metadata_path` must be NULL or a single non-empty character string.", call. = FALSE)
    }

    if (!file.exists(model_metadata_path)) {
      stop("Model metadata file does not exist: ", model_metadata_path, call. = FALSE)
    }
  }

  if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    stop("`assay_name` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit)) {
    stop("`unit` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.logical(normalized) || length(normalized) != 1L || is.na(normalized)) {
    stop("`normalized` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.character(id_col) || length(id_col) != 1L || is.na(id_col) || !nzchar(id_col)) {
    stop("`id_col` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(feature_type) || length(feature_type) != 1L || is.na(feature_type) || !nzchar(feature_type)) {
    stop("`feature_type` must be a single non-empty character string.", call. = FALSE)
  }

  ## ---------------------------------------------------
  ## Read expression file
  ## ---------------------------------------------------

  message("Reading protein expression assay file: ", basename(path))

  expression_df <- data.table::fread(
    file = path,
    check.names = FALSE,
    data.table = TRUE,
    showProgress = interactive()
  )

  if (nrow(expression_df) == 0L) {
    stop("No data was read from protein expression assay file.", call. = FALSE)
  }

  if (!id_col %in% names(expression_df)) {
    stop("Specified `id_col` not found in protein expression assay file: ", id_col, call. = FALSE)
  }

  message(
    "Successfully read protein expression assay file with ",
    nrow(expression_df),
    " rows and ",
    ncol(expression_df),
    " columns."
  )

  ## ---------------------------------------------------
  ## Filter to default model entries
  ## ---------------------------------------------------

  n_rows_before_default_filter <- nrow(expression_df)

  if ("IsDefaultEntryForModel" %in% names(expression_df)) {
    expression_df <- expression_df[
      IsDefaultEntryForModel == "Yes"
    ]
  } else {
    warning(
      "`IsDefaultEntryForModel` column was not found. ",
      "Protein expression loader could not filter to default model entries.",
      call. = FALSE
    )
  }

  n_rows_after_default_filter <- nrow(expression_df)

  if (n_rows_after_default_filter == 0L) {
    stop("No rows remained after filtering to default model entries.", call. = FALSE)
  }

  message(
    "Filtered protein expression assay from ",
    n_rows_before_default_filter,
    " to ",
    n_rows_after_default_filter,
    " default model entries."
  )

  ## ---------------------------------------------------
  ## Extract model identifiers
  ## ---------------------------------------------------

  model_ids <- expression_df[[id_col]]

  if (anyNA(model_ids) || any(!nzchar(model_ids))) {
    stop("`", id_col, "` contains NA or empty values.", call. = FALSE)
  }

  if (anyDuplicated(model_ids)) {
    duplicated_ids <- unique(model_ids[duplicated(model_ids)])

    stop(
      "`",
      id_col,
      "` contains duplicated values after filtering to default model entries. ",
      "Example duplicated IDs: ",
      paste(utils::head(duplicated_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Separate index columns from expression columns
  ## ---------------------------------------------------

  index_like_cols <- names(expression_df)[
    names(expression_df) == "" |
      names(expression_df) %in% c("V1")
  ]

  if (length(index_like_cols) > 0L) {
    message(
      "Dropping index-like column(s): ",
      paste(index_like_cols, collapse = ", ")
    )

    expression_df[
      ,
      (index_like_cols) := NULL
    ]
  }

  identifier_cols <- c(
    id_col,
    "SequencingID",
    "ModelConditionID",
    "ModelID",
    "IsDefaultEntryForMC",
    "IsDefaultEntryForModel"
  )

  identifier_cols <- intersect(identifier_cols, names(expression_df))
  expression_cols <- setdiff(names(expression_df), identifier_cols)

  if (length(expression_cols) == 0L) {
    stop(
      "No protein expression value columns were found after removing identifier columns.",
      call. = FALSE
    )
  }

  if (anyNA(expression_cols) || any(!nzchar(expression_cols))) {
    stop(
      "Protein expression feature columns contain NA or empty names after removing identifier columns.",
      call. = FALSE
    )
  }

  if (anyDuplicated(expression_cols)) {
    duplicated_features <- unique(expression_cols[duplicated(expression_cols)])

    stop(
      "Protein expression feature columns contain duplicated names. Examples: ",
      paste(utils::head(duplicated_features, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  expression_index <- as.data.frame(
    expression_df[
      ,
      ..identifier_cols
    ],
    stringsAsFactors = FALSE
  )

  rownames(expression_index) <- model_ids

  ## ---------------------------------------------------
  ## Validate expression columns
  ## ---------------------------------------------------

  if (length(identifier_cols) > 0L) {
    expression_df[
      ,
      (identifier_cols) := NULL
    ]
  }

  non_numeric_cols <- names(expression_df)[
    !vapply(expression_df, is.numeric, logical(1))
  ]

  if (length(non_numeric_cols) > 0L) {
    stop(
      "Some protein expression columns are not numeric. Examples: ",
      paste(utils::head(non_numeric_cols, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Convert expression values to matrix
  ## ---------------------------------------------------

  message(
    "Converting protein expression data to matrix in chunks: ",
    n_rows_after_default_filter,
    " models x ",
    length(expression_cols),
    " features."
  )

  expression_matrix <- matrix(
    NA_real_,
    nrow = length(expression_cols),
    ncol = n_rows_after_default_filter,
    dimnames = list(expression_cols, model_ids)
  )

  chunk_size <- 1000L
  chunk_starts <- seq.int(1L, length(expression_cols), by = chunk_size)

  for (chunk_start in chunk_starts) {
    chunk_end <- min(chunk_start + chunk_size - 1L, length(expression_cols))
    chunk_cols <- expression_cols[chunk_start:chunk_end]

    chunk_matrix <- as.matrix(
      expression_df[
        ,
        ..chunk_cols
      ]
    )

    expression_matrix[chunk_cols, ] <- t(chunk_matrix)
  }

  rm(expression_df)
  invisible(gc())

  if (!is.numeric(expression_matrix)) {
    stop("Protein expression matrix is not numeric after conversion.", call. = FALSE)
  }

  if (anyNA(expression_matrix)) {
    stop(
      "Protein expression matrix contains NA values after conversion.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Validate matrix dimnames
  ## ---------------------------------------------------

  feature_ids <- rownames(expression_matrix)
  model_ids <- colnames(expression_matrix)

  if (is.null(feature_ids) || anyNA(feature_ids) || any(!nzchar(feature_ids))) {
    stop(
      "Protein expression matrix feature IDs are missing, NA, or empty.",
      call. = FALSE
    )
  }

  if (anyDuplicated(feature_ids)) {
    duplicated_features <- unique(feature_ids[duplicated(feature_ids)])

    stop(
      "Protein expression matrix feature IDs contain duplicates. Examples: ",
      paste(utils::head(duplicated_features, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(model_ids) || anyNA(model_ids) || any(!nzchar(model_ids))) {
    stop(
      "Protein expression matrix model IDs are missing, NA, or empty.",
      call. = FALSE
    )
  }

  if (anyDuplicated(model_ids)) {
    duplicated_models <- unique(model_ids[duplicated(model_ids)])

    stop(
      "Protein expression matrix model IDs contain duplicates. Examples: ",
      paste(utils::head(duplicated_models, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Construct rowData and colData
  ## ---------------------------------------------------

  feature_metadata <- S4Vectors::DataFrame(
    feature_id = feature_ids,
    feature_type = feature_type,
    row.names = feature_ids
  )

  if (is.null(model_metadata_path)) {
    model_metadata_path <- .find_sidecar_file(
      path = path,
      filename = "Model.csv"
    )
  }

  model_metadata <- .load_model_metadata(
    model_metadata_path = model_metadata_path,
    model_ids = model_ids
  )

  ## Add the expression-file index columns to colData if they align.
  match_idx <- match(model_ids, rownames(expression_index))

  if (anyNA(match_idx)) {
    stop(
      "Could not align protein expression index metadata to expression matrix columns.",
      call. = FALSE
    )
  }

  expression_index <- expression_index[match_idx, , drop = FALSE]

  for (col in colnames(expression_index)) {
    if (!col %in% colnames(model_metadata)) {
      model_metadata[[col]] <- expression_index[[col]]
    }
  }

  ## ---------------------------------------------------
  ## Final consistency checks before ProteinExpressionAssay
  ## ---------------------------------------------------

  if (!identical(rownames(expression_matrix), rownames(feature_metadata))) {
    stop(
      "`rownames(expression_matrix)` and `rownames(feature_metadata)` do not match.",
      call. = FALSE
    )
  }

  if (!identical(colnames(expression_matrix), rownames(model_metadata))) {
    stop(
      "`colnames(expression_matrix)` and `rownames(model_metadata)` do not match.",
      call. = FALSE
    )
  }

  if (nrow(feature_metadata) != nrow(expression_matrix)) {
    stop(
      "`feature_metadata` must have one row per expression matrix row.",
      call. = FALSE
    )
  }

  if (nrow(model_metadata) != ncol(expression_matrix)) {
    stop(
      "`model_metadata` must have one row per expression matrix column.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Return ProteinExpressionAssay
  ## ---------------------------------------------------

  ProteinExpressionAssay(
    data = expression_matrix,
    rowData = feature_metadata,
    colData = model_metadata,
    metadata = list(
      assay_name = assay_name,
      source = "DepMap",
      source_file = normalizePath(path, winslash = "/", mustWork = TRUE),
      unit = unit,
      normalized = normalized,
      feature_type = feature_type,
      id_col = id_col,
      filtered_to_default_model_entries = "IsDefaultEntryForModel" %in% colnames(expression_index),
      n_rows_before_default_filter = n_rows_before_default_filter,
      n_rows_after_default_filter = n_rows_after_default_filter,
      model_metadata_loaded = !is.null(model_metadata_path),
      model_metadata_file = if (!is.null(model_metadata_path)) {
        normalizePath(model_metadata_path, winslash = "/", mustWork = TRUE)
      } else {
        NA_character_
      }
    ),
    assay_name = "protein_expression",
    unit = unit,
    normalized = normalized,
    feature_type = feature_type,
    source = "DepMap",
    source_file = normalizePath(path, winslash = "/", mustWork = TRUE)
  )
}
