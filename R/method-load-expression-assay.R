# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(data.table)

#' @keywords internal
#' @return An `ExpressionAssay`, which inherits from `OncoAssay`
load_expression_assay <- function(
  path,
  assay_name = "Expression",
  model_metadata_path = NULL,
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
    stop("Expression assay file does not exist: ", path, call. = FALSE)
  }

  if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    stop("`assay_name` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.null(model_metadata_path)) {
    if (!is.character(model_metadata_path) || length(model_metadata_path) != 1L || is.na(model_metadata_path) || !nzchar(model_metadata_path)) {
      stop("`model_metadata_path` must be NULL or a single non-empty character string.", call. = FALSE)
    }

    if (!file.exists(model_metadata_path)) {
      stop("Model metadata file does not exist: ", model_metadata_path, call. = FALSE)
    }
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
  header <- names(data.table::fread(
    file = path,
    nrows = 0L,
    check.names = FALSE,
    showProgress = FALSE
  ))


  identifier_cols <- c(
    "SequencingID",
    "ModelConditionID",
    "ModelID",
    "IsDefaultEntryForMC",
    "IsDefaultEntryForModel"
  )

  identifier_cols <- intersect(identifier_cols, header)

  if (!id_col %in% header) {
    stop("Specified `id_col` not found in expression assay file: ", id_col, call. = FALSE)
  }

  expression_cols <- setdiff(header, identifier_cols)

  if (length(expression_cols) == 0L) {
    stop("No expression columns were found in expression assay file.", call. = FALSE)
  }

  message(
    "Reading expression assay file with ",
    length(identifier_cols),
    " identifier columns and ",
    length(expression_cols),
    " expression columns."
  )

  expression_dt <- data.table::fread(
    file = path,
    check.names = FALSE,
    data.table = TRUE,
    showProgress = TRUE, # TODO: set to interactive() when done developing
    nThread = 2L,
    quote = "",
    na.strings = c("", "NA", "NaN", " "),
    colClasses = list(
      character = identifier_cols,
      numeric = expression_cols
    )
  )

  if (nrow(expression_dt) == 0L) {
    stop("No data was read from expression assay file.", call. = FALSE)
  }

  message(
    "Successfully read expression assay file with ",
    nrow(expression_dt),
    " rows and ",
    ncol(expression_dt),
    " columns."
  )
}
