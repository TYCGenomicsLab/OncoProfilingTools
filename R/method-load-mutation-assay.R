# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(data.table)

#' Load a DepMap mutation assay
#'
#' Loads the DepMap OmicsSomaticMutations.csv file into a MutationAssay.
#'
#' The mutation file is a MAF-like long table where each row is a somatic
#' mutation call. The loader keeps the call annotations in rowData, aligns model
#' metadata from Model.csv into colData, and stores the allele frequency values
#' in a model-by-call assay matrix.
#'
#' If Model.csv exists beside the mutation file, it is loaded automatically.
#'
#' @param path A character string specifying the path to the mutation assay file.
#' @param assay_name A character string specifying the assay name used in metadata.
#' @param model_metadata_path Optional path to Model.csv. If NULL, the loader
#'   searches for Model.csv beside path.
#' @param value_col A character string specifying the numeric source column to
#'   store in the assay matrix.
#'
#' @return A MutationAssay, which inherits from SummarizedExperiment.
#'
#' @include internal-metadata.R
#' @keywords internal
load_mutation_assay <- function(
  path,
  assay_name = "Mutation",
  model_metadata_path = NULL,
  value_col = "AF"
) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop("Mutation assay file does not exist: ", path, call. = FALSE)
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

  if (!is.character(value_col) || length(value_col) != 1L || is.na(value_col) || !nzchar(value_col)) {
    stop("`value_col` must be a single non-empty character string.", call. = FALSE)
  }

  message("Reading mutation assay file: ", basename(path))

  mutation_df <- data.table::fread(
    file = path,
    check.names = FALSE,
    data.table = TRUE,
    showProgress = interactive()
  )

  if (nrow(mutation_df) == 0L) {
    stop("No data was read from mutation assay file.", call. = FALSE)
  }

  if (!"ModelID" %in% names(mutation_df)) {
    stop("Mutation assay file must contain a `ModelID` column.", call. = FALSE)
  }

  if (!value_col %in% names(mutation_df)) {
    stop("Mutation assay file must contain a `", value_col, "` column.", call. = FALSE)
  }

  if (!is.numeric(mutation_df[[value_col]])) {
    stop("Mutation assay column `", value_col, "` must be numeric.", call. = FALSE)
  }

  message(
    "Successfully read mutation assay file with ",
    nrow(mutation_df),
    " rows and ",
    ncol(mutation_df),
    " columns."
  )

  if ("IsDefaultEntryForModel" %in% names(mutation_df)) {
    n_rows_before_default_filter <- nrow(mutation_df)

    mutation_df <- mutation_df[mutation_df$IsDefaultEntryForModel == "Yes", , drop = FALSE]

    if (nrow(mutation_df) == 0L) {
      stop("No rows remained after filtering to default model entries.", call. = FALSE)
    }

    message(
      "Filtered mutation assay from ",
      n_rows_before_default_filter,
      " to ",
      nrow(mutation_df),
      " default model entries."
    )
  } else {
    warning(
      "`IsDefaultEntryForModel` column was not found. Mutation loader could not filter to default model entries.",
      call. = FALSE
    )
  }

  if (is.null(model_metadata_path)) {
    model_metadata_path <- .find_sidecar_file(
      path = path,
      filename = "Model.csv"
    )
  }

  model_ids <- unique(as.character(mutation_df$ModelID))

  if (anyNA(model_ids) || any(!nzchar(model_ids))) {
    stop("`ModelID` contains NA or empty values.", call. = FALSE)
  }

  if (anyDuplicated(model_ids)) {
    duplicated_models <- unique(model_ids[duplicated(model_ids)])

    stop(
      "`ModelID` contains duplicated values. Examples: ",
      paste(utils::head(duplicated_models, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  model_metadata <- .load_model_metadata(
    model_metadata_path = model_metadata_path,
    model_ids = model_ids
  )

  n_models <- length(model_ids)
  row_ids <- paste0("mutation_", seq_len(nrow(mutation_df)))

  feature_cols <- intersect(
    c(
      "Chrom",
      "Pos",
      "Ref",
      "Alt",
      "HugoSymbol",
      "VariantInfo",
      "DNAChange",
      "ProteinChange",
      "ModelID",
      "SequencingID"
    ),
    names(mutation_df)
  )

  if (length(feature_cols) > 0L) {
    feature_labels <- vapply(
      seq_len(nrow(mutation_df)),
      function(i) {
        values <- vapply(
          feature_cols,
          function(col) as.character(mutation_df[[col]][i]),
          character(1)
        )

        label <- paste(values, collapse = "|")

        if (!nzchar(label)) {
          row_ids[[i]]
        } else {
          label
        }
      },
      character(1)
    )

    row_ids <- make.unique(feature_labels, sep = "_")
  }

  mutation_values <- mutation_df[[value_col]]
  model_index <- match(as.character(mutation_df$ModelID), model_ids)

  if (anyNA(model_index)) {
    missing_models <- unique(as.character(mutation_df$ModelID)[is.na(model_index)])

    stop(
      "Some mutation rows refer to ModelID values that were not retained in the metadata alignment. Examples: ",
      paste(utils::head(missing_models, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  mutation_matrix <- matrix(
    0,
    nrow = nrow(mutation_df),
    ncol = n_models,
    dimnames = list(row_ids, model_ids)
  )

  mutation_matrix[cbind(seq_len(nrow(mutation_df)), model_index)] <- mutation_values

  row_metadata <- S4Vectors::DataFrame(mutation_df)
  rownames(row_metadata) <- row_ids

  assay_metadata <- list(
    assay_name = assay_name,
    source = "DepMap",
    source_file = normalizePath(path, winslash = "/", mustWork = TRUE),
    value_col = value_col,
    filtered_to_default_model_entries = "IsDefaultEntryForModel" %in% names(mutation_df),
    model_metadata_loaded = !is.null(model_metadata_path),
    model_metadata_file = if (!is.null(model_metadata_path)) {
      normalizePath(model_metadata_path, winslash = "/", mustWork = TRUE)
    } else {
      NA_character_
    },
    n_rows = nrow(mutation_df),
    n_models = n_models
  )

  MutationAssay(
    data = mutation_matrix,
    rowData = row_metadata,
    colData = model_metadata,
    metadata = assay_metadata,
    assay_name = "mutation"
  )
}
