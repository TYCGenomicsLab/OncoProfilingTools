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
  assay_name = "PRISM Primary Repurposing",
  compound_metadata_path = NULL,
  cell_line_metadata_path = NULL,
  unit = "logfoldchange"
) {
  ## ---------------------------------------------------
  ## Validate arguments
  ## ---------------------------------------------------

  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop("Drug response assay file does not exist: ", path, call. = FALSE)
  }

  if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    stop("`assay_name` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit)) {
    stop("`unit` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.null(compound_metadata_path)) {
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
  }

  if (!is.null(cell_line_metadata_path)) {
    if (!is.character(cell_line_metadata_path) || length(cell_line_metadata_path) != 1L || is.na(cell_line_metadata_path) || !nzchar(cell_line_metadata_path)) {
      stop(
        "`cell_line_metadata_path` must be NULL or a single non-empty character string.",
        call. = FALSE
      )
    }

    if (!file.exists(cell_line_metadata_path)) {
      stop(
        "Cell line metadata file does not exist: ",
        cell_line_metadata_path,
        call. = FALSE
      )
    }
  }

  ## ---------------------------------------------------
  ## Resolve PRISM sidecar metadata files if not provided
  ## ---------------------------------------------------

  if (is.null(compound_metadata_path)) {
    compound_metadata_path <- .find_sidecar_file(
      path = path,
      filename = "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"
    )
  }

  if (is.null(cell_line_metadata_path)) {
    cell_line_metadata_path <- .find_sidecar_file(
      path = path,
      filename = "Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv"
    )
  }

  ## ---------------------------------------------------
  ## Read response matrix
  ## ---------------------------------------------------

  message("Reading drug response assay file: ", basename(path))

  response_df <- data.table::fread(
    file = path,
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

  treatment_ids <- response_df[[1L]]
  model_ids <- colnames(response_df)[-1L]

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
    dplyr::select(-1L)

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

  treatment_metadata <- .load_prism_compound_metadata(
    compound_metadata_path = compound_metadata_path,
    treatment_ids = treatment_ids
  )

  model_metadata <- .load_prism_cell_line_metadata(
    cell_line_metadata_path = cell_line_metadata_path,
    model_ids = model_ids
  )

  ## ---------------------------------------------------
  ## Final consistency checks
  ## ---------------------------------------------------

  if (!identical(rownames(response_matrix), rownames(treatment_metadata))) {
    stop(
      "`rownames(response_matrix)` and `rownames(treatment_metadata)` do not match.",
      call. = FALSE
    )
  }

  if (!identical(colnames(response_matrix), rownames(model_metadata))) {
    stop(
      "`colnames(response_matrix)` and `rownames(model_metadata)` do not match.",
      call. = FALSE
    )
  }

  if (nrow(treatment_metadata) != nrow(response_matrix)) {
    stop(
      "`treatment_metadata` must have one row per response matrix row.",
      call. = FALSE
    )
  }

  if (nrow(model_metadata) != ncol(response_matrix)) {
    stop(
      "`model_metadata` must have one row per response matrix column.",
      call. = FALSE
    )
  }

  ## ---------------------------------------------------
  ## Return DrugResponseAssay
  ## ---------------------------------------------------

  assay_metadata <- list(
    assay_name = assay_name,
    source = "PRISM",
    source_file = normalizePath(path, winslash = "/", mustWork = TRUE),
    compound_metadata_file = if (is.null(compound_metadata_path)) {
      NA_character_
    } else {
      normalizePath(compound_metadata_path, winslash = "/", mustWork = TRUE)
    },
    cell_line_metadata_file = if (is.null(cell_line_metadata_path)) {
      NA_character_
    } else {
      normalizePath(cell_line_metadata_path, winslash = "/", mustWork = TRUE)
    },
    unit = unit,
    n_treatments = nrow(response_matrix),
    n_models = ncol(response_matrix),
    n_missing_values = n_missing,
    missing_value_percent = if (length(response_matrix) == 0L) {
      NA_real_
    } else {
      round(100 * n_missing / length(response_matrix), digits = 2L)
    },
    compound_metadata_loaded = !is.null(compound_metadata_path),
    cell_line_metadata_loaded = !is.null(cell_line_metadata_path)
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(
      response = response_matrix
    ),
    rowData = treatment_metadata,
    colData = model_metadata,
    metadata = assay_metadata
  )

  methods::as(se, "DrugResponseAssay")
}
