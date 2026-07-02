# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(dplyr)
library(enumr)

utils::globalVariables(c(
  "mean_diff",
  "neg_log10_p",
  "significant",
  "feature_label"
))

#' Drug sensitivity comparison result
#'
#' @keywords internal
methods::setClass(
  "DrugSensitivityComparison",
  slots = c(
    stats = "data.frame",
    plot_data = "list",
    group1_n_models = "integer",
    group2_n_models = "integer",
    unit = "character",
    p_adj_method = "character",
    effect_threshold = "numeric",
    group1_label = "character",
    group2_label = "character"
  )
)

#' Drug sensitivity response units
#' LFC = log fold change, 
#' AUC = area under the curve,
#' IC50 = half maximal inhibitory concentration
#' @export
SensitivityUnits <- enumr::enum(
  "LFC" = "logfoldchange",
  # TODO: "AUC" = "AUC",
  # TODO: "IC50" = "IC50"
)



#' Compare drug sensitivities between two OncoExperiment groups
#'
#' Compares PRISM drug response values between two user-defined groups stored in
#' a single `OncoExperiment` object using per-drug Welch two-sample t-tests.
#' 
#' **User is expected to define the groups prior to calling this method.**
#' For example, to compare drug sensitivities between two cancer types, subset
#' the `OncoExperiment` and then assign groups:
#' `COAD$group <- c(1, 1, 2, 1, 2, 2, 1, ...)`
#' 
#' Whatever labels were chosen for the groups must be passed into the function.
#' The method looks in the top-level `$group` column and compares the selected
#' values within the PRISM assay.
#'
#' This method expects a loadable `DrugResponseAssay` (typically named
#' "PRISM") with assay data named "response".
#'
#' The sign convention is:
#' `mean_diff = mean(group1) - mean(group2)`.
#' Negative values indicate lower average response in `group1`.
#'
#' @param object An `OncoExperiment` object with a top-level `group` column.
#' @param group1 Label for the first group of models/samples.
#' @param group2 Label for the second group of models/samples.
#' @param unit Optional character scalar describing the expected response unit
#'   (for example, `"LFC"`). If `NULL`, the unit is inferred from the assay
#'   metadata when available and normalized to `"LFC"` when the metadata uses
#'   an equivalent label such as `"logfoldchange"`.
#' @param p_adj_method Method passed to `stats::p.adjust()`.
#' @param effect_threshold Non-negative numeric threshold for
#'   `abs(mean_diff)` when computing significance flags.
#' @param ... Additional arguments. Currently unused.
#'
#' @return A `DrugSensitivityComparison` object.
#' @export
methods::setMethod(
  f = "compare_drug_sensitivities",
  signature = c(object = "OncoExperiment", group1 = "ANY", group2 = "ANY"),
  definition = function(object,
                        group1,
                        group2,
                        unit = NULL,
                        p_adj_method = "BH",
                        effect_threshold = 0,
                        ...) {
    if (!is.atomic(group1) || length(group1) != 1L || is.na(group1)) {
      stop("`group1` must be a single non-NA value.", call. = FALSE)
    }

    if (!is.atomic(group2) || length(group2) != 1L || is.na(group2)) {
      stop("`group2` must be a single non-NA value.", call. = FALSE)
    }

    group1_key <- as.character(group1)
    group2_key <- as.character(group2)

    if (!nzchar(group1_key) || !nzchar(group2_key)) {
      stop("`group1` and `group2` must be non-empty labels.", call. = FALSE)
    }

    if (identical(group1_key, group2_key)) {
      stop("`group1` and `group2` must be different labels.", call. = FALSE)
    }

    if (!is.null(unit) && (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit))) {
      stop("`unit` must be a single non-empty character string.", call. = FALSE)
    }

    if (!is.character(p_adj_method) || length(p_adj_method) != 1L || is.na(p_adj_method) || !nzchar(p_adj_method)) {
      stop("`p_adj_method` must be a single non-empty character string.", call. = FALSE)
    }

    if (!is.numeric(effect_threshold) || length(effect_threshold) != 1L || is.na(effect_threshold) || effect_threshold < 0) {
      stop("`effect_threshold` must be a single non-negative numeric value.", call. = FALSE)
    }

    .normalize_response_unit <- function(x) {
      if (is.null(x)) {
        return(NULL)
      }

      x_chr <- tolower(as.character(x))

      if (x_chr %in% c("logfoldchange", "log fold change", "lfc")) {
        return("LFC")
      }

      as.character(x)
    }

    .get_drug_response_assay <- function(x) {
      experiments <- MultiAssayExperiment::experiments(x)

      if ("PRISM" %in% names(experiments) && methods::is(experiments[["PRISM"]], "DrugResponseAssay")) {
        return(list(assay = experiments[["PRISM"]], assay_name = "PRISM"))
      }

      drug_assay_idx <- vapply(
        experiments,
        function(experiment) methods::is(experiment, "DrugResponseAssay"),
        logical(1)
      )

      if (!any(drug_assay_idx)) {
        stop(
          "No `DrugResponseAssay` found in the provided OncoExperiment. ",
          "Load PRISM first with `load_assays(..., assay_type = \"PRISM\")`.",
          call. = FALSE
        )
      }

      assay_idx <- which(drug_assay_idx)[1L]
      list(
        assay = experiments[[assay_idx]],
        assay_name = names(experiments)[assay_idx]
      )
    }

    assay_info <- .get_drug_response_assay(object)
    assay <- assay_info$assay
    assay_name <- assay_info$assay_name

    assay_unit <- .normalize_response_unit(S4Vectors::metadata(assay)$unit)

    unit <- .normalize_response_unit(unit)

    if (!is.null(unit) && !is.null(assay_unit) && !identical(assay_unit, unit)) {
      stop(
        "`unit` does not match the loaded assay metadata unit. Expected `",
        assay_unit,
        "`, got `",
        unit,
        "`.",
        call. = FALSE
      )
    }

    resolved_unit <- if (!is.null(assay_unit)) assay_unit else unit

    if (is.null(resolved_unit)) {
      stop(
        "Unable to resolve the response unit. Pass `unit` or load assay metadata with a unit.",
        call. = FALSE
      )
    }

    response <- SummarizedExperiment::assay(
      assay,
      i = "response",
      withDimnames = TRUE
    )

    if (!is.numeric(response)) {
      stop("The response assay must be a numeric matrix.", call. = FALSE)
    }

    col_data <- MultiAssayExperiment::colData(object)

    if (!"group" %in% colnames(col_data)) {
      stop("`object` must contain a `group` column in `colData(object)`.", call. = FALSE)
    }

    if (nrow(col_data) == 0L) {
      stop("`object` has no models/samples in `colData(object)`.", call. = FALSE)
    }

    sample_map <- MultiAssayExperiment::sampleMap(object)
    assay_map <- sample_map[sample_map$assay == assay_name, , drop = FALSE]

    if (nrow(assay_map) == 0L) {
      stop(
        "No sample-map entries were found for the loaded drug response assay. ",
        "Ensure the assay is aligned to the top-level `group` column.",
        call. = FALSE
      )
    }

    primary_lookup <- as.character(MultiAssayExperiment::colData(object)$group)
    names(primary_lookup) <- rownames(col_data)

    if (is.null(rownames(col_data)) || anyNA(rownames(col_data))) {
      stop("`colData(object)` must have valid row names for group matching.", call. = FALSE)
    }

    assay_primary_ids <- assay_map$primary[match(colnames(response), as.character(assay_map$colname))]

    if (anyNA(assay_primary_ids)) {
      stop(
        "Could not map all drug response assay columns back to the top-level `colData(object)`.",
        call. = FALSE
      )
    }

    assay_groups <- primary_lookup[assay_primary_ids]

    assigned_idx <- !is.na(assay_groups) & nzchar(assay_groups)

    if (any(!assigned_idx)) {
      message(
        "Ignoring ",
        sum(!assigned_idx),
        " assay sample(s) without a group assignment in `colData(object)`."
      )
      assay_groups <- assay_groups[assigned_idx]
      response <- response[, assigned_idx, drop = FALSE]
    }

    if (ncol(response) == 0L) {
      stop(
        "No assay samples remain after removing unassigned group labels.",
        call. = FALSE
      )
    }

    group1_idx <- which(assay_groups == group1_key)
    group2_idx <- which(assay_groups == group2_key)

    if (length(group1_idx) == 0L || length(group2_idx) == 0L) {
      stop(
        "The requested group labels were not found in the top-level `group` column.",
        call. = FALSE
      )
    }

    response1 <- response[, group1_idx, drop = FALSE]
    response2 <- response[, group2_idx, drop = FALSE]

    if (ncol(response1) == 0L || ncol(response2) == 0L) {
      stop(
        "Both groups must contain at least one model/sample in the PRISM assay.",
        call. = FALSE
      )
    }

    common_features <- rownames(response)

    n_group1 <- rowSums(!is.na(response1))
    n_group2 <- rowSums(!is.na(response2))

    mean_group1 <- rowMeans(response1, na.rm = TRUE)
    mean_group2 <- rowMeans(response2, na.rm = TRUE)

    mean_group1[is.nan(mean_group1)] <- NA_real_
    mean_group2[is.nan(mean_group2)] <- NA_real_

    mean_diff <- mean_group1 - mean_group2

    p_values <- vapply(
      X = seq_len(nrow(response1)),
      FUN = function(i) {
        x <- response1[i, ]
        y <- response2[i, ]

        x <- x[is.finite(x)]
        y <- y[is.finite(y)]

        if (length(x) < 2L || length(y) < 2L) {
          return(NA_real_)
        }

        tryCatch(
          stats::t.test(x, y, var.equal = FALSE)$p.value,
          error = function(e) NA_real_
        )
      },
      FUN.VALUE = numeric(1)
    )

    p_adj <- stats::p.adjust(p_values, method = p_adj_method)

    results <- data.frame(
      feature_id = common_features,
      n_group1 = as.integer(n_group1),
      n_group2 = as.integer(n_group2),
      mean_group1 = as.numeric(mean_group1),
      mean_group2 = as.numeric(mean_group2),
      mean_diff = as.numeric(mean_diff),
      p_value = as.numeric(p_values),
      p_adj = as.numeric(p_adj),
      stringsAsFactors = FALSE
    )

    row_metadata <- as.data.frame(
      SummarizedExperiment::rowData(assay),
      stringsAsFactors = FALSE
    )

    if (nrow(row_metadata) > 0L) {
      row_metadata$feature_id <- rownames(row_metadata)
      row_metadata <- row_metadata[row_metadata$feature_id %in% common_features, , drop = FALSE]

      if (nrow(row_metadata) > 0L) {
        keep_cols <- c(
          "feature_id",
          "treatment_id",
          "compound_id",
          "drug_name",
          "moa",
          "target",
          "IDs",
          "Drug.Name",
          "MOA",
          "repurposing_target"
        )

        keep_cols <- intersect(keep_cols, colnames(row_metadata))

        if (length(keep_cols) > 1L) {
          row_metadata <- row_metadata[, keep_cols, drop = FALSE]
          results <- dplyr::left_join(results, row_metadata, by = "feature_id")
        }
      }
    }

    results <- results |>
      dplyr::mutate(
        neg_log10_p = -log10(.data$p_value),
        significant = !is.na(.data$p_adj) & .data$p_adj < 0.05 & abs(.data$mean_diff) > effect_threshold
      ) |>
      dplyr::arrange(.data$p_value)

    methods::new(
      "DrugSensitivityComparison",
      stats = results,
      plot_data = list(
        response = response,
        sample_groups = assay_groups,
        feature_ids = common_features,
        group1_idx = group1_idx,
        group2_idx = group2_idx,
        assay_name = assay_name,
        group1_label = group1_key,
        group2_label = group2_key
      ),
      group1_n_models = as.integer(ncol(response1)),
      group2_n_models = as.integer(ncol(response2)),
      unit = as.character(resolved_unit),
      p_adj_method = as.character(p_adj_method),
      effect_threshold = as.numeric(effect_threshold),
      group1_label = group1_key,
      group2_label = group2_key
    )
  }
)

#' Export a `DrugSensitivityComparison` object
#'
#' @param object A `DrugSensitivityComparison` object.
#' @param output_file Character scalar path to the `.xlsx` file to create.
#' @param overwrite Logical value; overwrite an existing file if `TRUE`.
#' @param include_volcano Logical value; if `TRUE`, add a volcano plot sheet.
#' @param n_label_volcano Integer scalar controlling how many features to label.
#' @param group1_label Optional override for the workbook label of group 1.
#' @param group2_label Optional override for the workbook label of group 2.
#' @param ... Additional arguments passed through to the workbook exporter.
#' @return Invisibly returns `output_file`.
#' @export
methods::setMethod(
  f = "export",
  signature = c(object = "DrugSensitivityComparison"),
  definition = function(
    object,
    output_file = "drug_sensitivities_comparison.xlsx",
    overwrite = TRUE,
    include_volcano = TRUE,
    n_label_volcano = 10L,
    group1_label = object@group1_label %||% "Group 1",
    group2_label = object@group2_label %||% "Group 2",
    ...
  ) {
    export_compare_drug_sensitivities_workbook(
      comparison_result = object,
      group1_label = group1_label,
      group2_label = group2_label,
      output_file = output_file,
      overwrite = overwrite,
      include_volcano = include_volcano,
      n_label_volcano = n_label_volcano
    )
  }
)

`%||%` <- function(a, b) if (!is.null(a)) a else b
