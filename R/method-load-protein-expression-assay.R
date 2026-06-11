# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(data.table)

#' Load a DepMap protein expression assay
#'
#' Loads a DepMap expression file into an `ProteinExpressionAssay`.
#'
#' The expression matrix is filtered to default model entries and stored as a
#' numeric matrix with `ModelID` values as row names and expression features as
#' columns. Metadata from `Model.csv` is not loaded here.
#'
#' @param path A character string specifying the path to the expression assay file.
#' @param assay_name A character string specifying the assay name.
#' @param unit A character string specifying the expression unit.
#' @param normalized A logical value indicating whether the expression data are normalized.
#' @param id_col A character string specifying the expression file column to use
#'   as row names.
#' @param feature_type A character string describing the expression feature type.
#'
#' @return An `ExpressionAssay`, which inherits from `OncoAssay`.
#' @keywords internal
load_protein_expression_assay <- function(
  path,
  assay_name = "Protein Expression",
  unit = "log2(TPM+1)",
  normalized = TRUE,
  id_col = "ModelID",
  feature_type = "transcript"
) {
  ## ---------------------------------------------------
  ## Validate arguments
  ## ---------------------------------------------------

  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop("Protein Expression assay file does not exist: ", path, call. = FALSE)
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

  message("Reading expression assay file: ", basename(path))

  expression_df <- data.table::fread(
    file = path,
    check.names = FALSE,
    data.table = TRUE,
    showProgress = interactive()
  )

  if (nrow(expression_df) == 0L) {
    stop("No data was read from expression assay file.", call. = FALSE)
  }

  if (!id_col %in% names(expression_df)) {
    stop("Specified `id_col` not found in expression assay file: ", id_col, call. = FALSE)
  }

  message(
    "Successfully read expression assay file with ",
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
    expression_df <- dplyr::filter(
      expression_df,
      .data$IsDefaultEntryForModel == "Yes"
    )
  } else {
    warning(
      "`IsDefaultEntryForModel` column was not found. ",
      "Expression loader could not filter to default model entries.",
      call. = FALSE
    )
  }

  n_rows_after_default_filter <- nrow(expression_df)

  if (n_rows_after_default_filter == 0L) {
    stop("No rows remained after filtering to default model entries.", call. = FALSE)
  }

  message(
    "Filtered expression assay from ",
    n_rows_before_default_filter,
    " to ",
    n_rows_after_default_filter,
    " default model entries."
  )

  ## ---------------------------------------------------
  ## Extract model identifiers
  ## ---------------------------------------------------

  row_ids <- expression_df[[id_col]]

  if (anyNA(row_ids) || any(!nzchar(row_ids))) {
    stop("`", id_col, "` contains NA or empty values.", call. = FALSE)
  }

  if (anyDuplicated(row_ids)) {
    duplicated_ids <- unique(row_ids[duplicated(row_ids)])

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

  identifier_cols <- c(
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
      "No expression value columns were found after removing identifier columns.",
      call. = FALSE
    )
  }

  expression_index <- expression_df |>
    dplyr::select(dplyr::all_of(identifier_cols)) |>
    as.data.frame(stringsAsFactors = FALSE)

  rownames(expression_index) <- row_ids

  ## ---------------------------------------------------
  ## Validate expression columns
  ## ---------------------------------------------------

  expression_values_df <- expression_df |>
    dplyr::select(dplyr::all_of(expression_cols))

  non_numeric_cols <- names(expression_values_df)[
    !vapply(expression_values_df, is.numeric, logical(1))
  ]

  if (length(non_numeric_cols) > 0L) {
    stop(
      "Some expression columns are not numeric. Examples: ",
      paste(utils::head(non_numeric_cols, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Convert expression values to matrix
  ## ---------------------------------------------------

  message(
    "Converting expression data to matrix: ",
    n_rows_after_default_filter,
    " rows x ",
    length(expression_cols),
    " features."
  )

  expression_matrix <- as.matrix(expression_values_df)
  rownames(expression_matrix) <- row_ids

  rm(expression_df, expression_values_df)
  invisible(gc())

  if (!is.numeric(expression_matrix)) {
    stop("Expression matrix is not numeric after conversion.", call. = FALSE)
  }

  if (anyNA(expression_matrix)) {
    stop(
      "Expression matrix contains NA values after conversion.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Construct minimal row and column metadata
  ## ---------------------------------------------------

  model_metadata <- data.frame(
    ModelID = row_ids,
    stringsAsFactors = FALSE
  )
  rownames(model_metadata) <- row_ids

  feature_ids <- colnames(expression_matrix)

  feature_metadata <- data.frame(
    feature_id = feature_ids,
    feature_type = feature_type,
    stringsAsFactors = FALSE
  )
  rownames(feature_metadata) <- feature_ids

  ## ---------------------------------------------------
  ## Return ExpressionAssay
  ## ---------------------------------------------------

  ProteinExpressionAssay(
    data = expression_matrix,
    metadata = list(
      source = "DepMap",
      source_file = normalizePath(path, winslash = "/", mustWork = TRUE),
      unit = unit,
      normalized = normalized,
      feature_type = feature_type,
      id_col = id_col,
      filtered_to_default_model_entries = "IsDefaultEntryForModel" %in% colnames(expression_index),
      n_rows_before_default_filter = n_rows_before_default_filter,
      n_rows_after_default_filter = n_rows_after_default_filter,
      model_metadata_loaded = FALSE
    ),
    row_data = list(
      models = model_metadata,
      expression_index = expression_index
    ),
    col_data = list(
      features = feature_metadata
    ),
    layers = list(),
    name = assay_name,
    unit = unit,
    normalized = normalized
  )
}