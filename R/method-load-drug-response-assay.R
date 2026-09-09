# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(data.table)

#' Load a PRISM drug response assay
#'
#' Loads a PRISM primary repurposing response matrix into a `DrugResponseAssay`.
#'
#' The PRISM response matrix is expected to contain treatments/compounds as rows
#' and DepMap model IDs as columns. This already matches the Bioconductor
#' convention used by `SummarizedExperiment`: rows are features and columns are
#' samples/models.
#'
#' Compound metadata from `Extended_Primary_Compound_List.csv` is stored in
#' `rowData()`. Cell line metadata from `Cell_Line_Meta_Data.csv` is stored in
#' `colData()`.
#'
#' @param path A character string specifying the file path to the drug response
#'   assay data.
#' @param assay_name A character string specifying the assay name.
#' @param compound_metadata_path Optional path to PRISM compound metadata.
#' @param cell_line_metadata_path Optional path to PRISM cell line metadata.
#' @param unit A character string specifying the unit of measurement.
#'
#' @return A `DrugResponseAssay`, which inherits from `SummarizedExperiment`.
#'
#' @include internal-metadata.R
#' @keywords internal
load_drug_response_assay <- function(
  data_path,
  metadata_path = NULL,
  assay_name = "DrugResponse",
  unit = "logfoldchange",
  identifier = NULL,
  ...
) {
  ## ---------------------------------------------------
  ## Validate arguments
  ## ---------------------------------------------------

  if (!is.character(data_path) || length(data_path) != 1L || is.na(data_path) || !nzchar(data_path)) {
    stop("`data_path` must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(data_path)) {
    stop("Drug response assay file does not exist: ", data_path, call. = FALSE)
  }

  if (!is.null(metadata_path)) {
    if (!is.character(metadata_path) || length(metadata_path) != 1L || is.na(metadata_path) || !nzchar(metadata_path)) {
      stop("`metadata_path` must be NULL or a single non-empty character string.", call. = FALSE)
    }

    if (!file.exists(metadata_path)) {
      stop("Model metadata file does not exist: ", metadata_path, call. = FALSE)
    }
  }

  if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    stop("`assay_name` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit)) {
    stop("`unit` must be a single non-empty character string.", call. = FALSE)
  }

  ## ---------------------------------------------------
  ## Read response matrix
  ## ---------------------------------------------------

  message("Reading drug response assay file: ", basename(data_path))

  response_df <- data.table::fread(
    file = data_path,
    check.names = FALSE,
    data.table = TRUE,
    showProgress = interactive()
  )

  if (nrow(response_df) == 0L) {
    stop("No data was read from drug response assay file.", call. = FALSE)
  }

  if (ncol(response_df) < 2L) {
    stop(
      "Drug response assay file must contain one treatment ID column and at least one model column.",
      call. = FALSE
    )
  }

  message(
    "Successfully read drug response assay file with ",
    nrow(response_df),
    " rows and ",
    ncol(response_df),
    " columns."
  )

  ## ---------------------------------------------------
  ## Extract treatment and model IDs
  ## ---------------------------------------------------

  # if user supplied identifier, that should be the colum name of ID's
  # else we assume it is in the first column
  if (!is.null(identifier)) {
    if (!is.character(identifier) || length(identifier) != 1L || is.na(identifier) || !nzchar(identifier)) {
      stop("`identifier` must be a single non-empty character string.", call. = FALSE)
    }
    if (!identifier %in% colnames(response_df)) {
      stop("Identifier column '", identifier, "' not found in data file: ", data_path, call. = FALSE)
    }
    treatment_ids <- response_df[[identifier]]
  } else {
    treatment_ids <- response_df[[1L]]
  }

  id_col <- if (is.null(identifier)) colnames(response_df)[1L] else identifier
  model_ids <- setdiff(colnames(response_df), id_col)

  if (anyNA(treatment_ids) || any(!nzchar(treatment_ids))) {
    stop(
      "Drug response treatment IDs contain NA or empty values.",
      call. = FALSE
    )
  }

  if (anyDuplicated(treatment_ids)) {
    duplicated_treatments <- unique(treatment_ids[duplicated(treatment_ids)])

    stop(
      "Drug response treatment IDs contain duplicated values. Examples: ",
      paste(utils::head(duplicated_treatments, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(model_ids) || any(!nzchar(model_ids))) {
    stop(
      "Drug response model columns contain NA or empty names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(model_ids)) {
    duplicated_models <- unique(model_ids[duplicated(model_ids)])

    stop(
      "Drug response model columns contain duplicated names. Examples: ",
      paste(utils::head(duplicated_models, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Extract numeric response values
  ## ---------------------------------------------------

  response_values_df <- response_df |>
    dplyr::select(-dplyr::all_of(id_col))

  non_numeric_cols <- names(response_values_df)[
    !vapply(response_values_df, is.numeric, logical(1))
  ]

  if (length(non_numeric_cols) > 0L) {
    stop(
      "Some drug response columns are not numeric. Examples: ",
      paste(utils::head(non_numeric_cols, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Convert to matrix
  ## ---------------------------------------------------

  message(
    "Converting drug response data to matrix: ",
    length(treatment_ids),
    " treatments x ",
    length(model_ids),
    " models."
  )

  response_matrix <- as.matrix(response_values_df)

  rownames(response_matrix) <- treatment_ids
  colnames(response_matrix) <- model_ids

  rm(response_df, response_values_df)
  invisible(gc())

  if (!is.numeric(response_matrix)) {
    stop(
      "Drug response matrix is not numeric after conversion.",
      call. = FALSE
    )
  }

  n_missing <- sum(is.na(response_matrix))

  if (n_missing > 0L) {
    missing_percent <- round(
      100 * n_missing / length(response_matrix),
      digits = 2L
    )

    message(
      "Drug response matrix contains ",
      n_missing,
      " missing values (",
      missing_percent,
      "%)."
    )
  }

  ## ---------------------------------------------------
  ## DO NOT TRANSPOSE
  ## ---------------------------------------------------

  ## PRISM response matrix is already treatments x models.
  ## That is already the SummarizedExperiment convention:
  ## rows = features/treatments, columns = samples/models.

  treatment_ids <- rownames(response_matrix)
  model_ids <- colnames(response_matrix)

  ## ---------------------------------------------------
  ## Construct rowData and colData
  ## ---------------------------------------------------

  # treatment_metadata <- .load_prism_compound_metadata(
  #  compound_metadata_path = compound_metadata_path,
  #  treatment_ids = treatment_ids
  # )

  # model_metadata <- .load_prism_cell_line_metadata(
  #  cell_line_metadata_path = cell_line_metadata_path,
  #  model_ids = model_ids
  # )

  # TODO: update to take generic drug response metadata file using same identifier info
  # load in the metadata file
  if (is.null(metadata_path)) {
    rowData <- S4Vectors::DataFrame(
      treatment_id = treatment_ids,
      row.names = treatment_ids
    )
  } else {
    metadata_df <- data.table::fread(
      file = metadata_path,
      check.names = FALSE,
      data.table = TRUE,
      showProgress = interactive()
    )

    if (nrow(metadata_df) == 0L) {
      stop("Metadata file is empty: ", metadata_path, call. = FALSE)
    }

    if (!id_col %in% colnames(metadata_df)) {
      stop("Identifier column '", id_col, "' not found in metadata file: ", metadata_path, call. = FALSE)
    }

    metadata_ids <- as.character(metadata_df[[id_col]])

    if (anyNA(metadata_ids) || any(!nzchar(metadata_ids))) {
      stop(
        "Metadata identifiers contain NA or empty values in column '",
        id_col,
        "'.",
        call. = FALSE
      )
    }

    if (anyDuplicated(metadata_ids)) {
      stop("Metadata file contains duplicated identifiers in column '", id_col, "'.", call. = FALSE)
    }

    match_idx <- match(treatment_ids, metadata_ids)
    if (anyNA(match_idx)) {
      warning("Some treatment IDs in the response matrix do not have matching metadata entries.", call. = FALSE)
    }

    row_data_df <- metadata_df[match_idx, , drop = FALSE]
    rowData <- S4Vectors::DataFrame(row_data_df)
    rownames(rowData) <- treatment_ids
  }

  colData <- S4Vectors::DataFrame(
    model_id = model_ids,
    row.names = model_ids
  )

  ## ---------------------------------------------------
  ## Final consistency checks
  ## ---------------------------------------------------


  ## ---------------------------------------------------
  ## Return DrugResponseAssay
  ## ---------------------------------------------------

  assay_metadata <- list(
    assay_name = assay_name,
    source = "DrugResponse",
    source_file = normalizePath(data_path, winslash = "/", mustWork = TRUE),
    unit = unit,
    n_treatments = nrow(response_matrix),
    n_models = ncol(response_matrix),
    n_missing_values = n_missing,
    missing_value_percent = if (length(response_matrix) == 0L) {
      NA_real_
    } else {
      round(100 * n_missing / length(response_matrix), digits = 2L)
    }
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(
      response = response_matrix
    ),
    rowData = rowData,
    colData = colData,
    metadata = assay_metadata
  )

  methods::as(se, "DrugResponseAssay")
}
