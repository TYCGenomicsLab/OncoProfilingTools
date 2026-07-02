# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Compare drug sensitivities between two groups
#'
#' S4 generic for comparing PRISM drug sensitivities between two groups stored
#' in a single `OncoExperiment` object.
#'
#' @param object An `OncoExperiment` object containing a top-level `group`
#'   column.
#' @param group1 Label identifying group 1.
#' @param group2 Label identifying group 2.
#' @param unit Character scalar describing the expected response unit.
#' @param p_adj_method Character scalar passed to `stats::p.adjust()`.
#' @param effect_threshold Non-negative numeric effect-size cutoff.
#' @param ... Additional method-specific arguments.
#'
#' @return Method-specific comparison output.
#' @export
methods::setGeneric(
  name = "compare_drug_sensitivities",
  def = function(object, group1, group2, unit = NULL, p_adj_method = "BH", effect_threshold = 0, ...) {
    standardGeneric("compare_drug_sensitivities")
  },
  signature = c(
    "object",
    "group1",
    "group2",
    "unit",
    "p_adj_method",
    "effect_threshold"
  )
)

#' Export a drug sensitivity comparison to an Excel workbook
#'
#' Internal workbook writer used by the `export()` S4 method for
#' `DrugSensitivityComparison` objects.
#'
#' @param comparison_result A `DrugSensitivityComparison` object.
#' @param group1_label Character scalar used for display names in the workbook.
#' @param group2_label Character scalar used for display names in the workbook.
#' @param output_file Character scalar path to the `.xlsx` file to create.
#' @param overwrite Logical value; overwrite an existing file if `TRUE`.
#' @param include_volcano Logical value; if `TRUE`, add a volcano plot sheet.
#' @param n_label_volcano Integer scalar controlling how many features to label.
#' @return Invisibly returns `output_file`.
#' @keywords internal
export_compare_drug_sensitivities_workbook <- function(
  comparison_result,
  group1_label,
  group2_label,
  output_file = "drug_sensitivities_comparison.xlsx",
  overwrite = TRUE,
  include_volcano = TRUE,
  n_label_volcano = 10L
) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package `openxlsx` is required for workbook export.", call. = FALSE)
  }

  if (!methods::is(comparison_result, "DrugSensitivityComparison")) {
    stop("`comparison_result` must be a `DrugSensitivityComparison` object.", call. = FALSE)
  }

  if (!is.character(group1_label) || length(group1_label) != 1L || is.na(group1_label) || !nzchar(group1_label)) {
    stop("`group1_label` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(group2_label) || length(group2_label) != 1L || is.na(group2_label) || !nzchar(group2_label)) {
    stop("`group2_label` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(output_file) || length(output_file) != 1L || is.na(output_file) || !nzchar(output_file)) {
    stop("`output_file` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(include_volcano) || length(include_volcano) != 1L || is.na(include_volcano)) {
    stop("`include_volcano` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(n_label_volcano) || length(n_label_volcano) != 1L || is.na(n_label_volcano) || n_label_volcano < 0) {
    stop("`n_label_volcano` must be a non-negative numeric value.", call. = FALSE)
  }

  .clean_label <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) {
      x <- "Comparison"
    }
    x
  }

  .make_unique_sheet_name <- function(base_name, existing_names) {
    base_name <- .clean_label(base_name)
    base_name <- substr(base_name, 1L, 31L)
    candidate <- base_name
    suffix_idx <- 2L

    while (candidate %in% existing_names) {
      suffix <- paste0("_", suffix_idx)
      candidate <- paste0(substr(base_name, 1L, max(1L, 31L - nchar(suffix))), suffix)
      suffix_idx <- suffix_idx + 1L
    }

    candidate
  }

  .display_label <- function(x) {
    x <- as.character(x)
    ifelse(nzchar(x), x, "Unknown")
  }

  .feature_title <- function(row) {
    label <- row[["feature_label"]]
    if (!is.null(label) && !is.na(label) && nzchar(as.character(label))) {
      return(as.character(label))
    }
    as.character(row[["feature_id"]])
  }

  .make_plot_sheet <- function(wb, sheet_name, plot_fun, width = 8.5, height = 6.5) {
    openxlsx::addWorksheet(wb, sheet_name, gridLines = FALSE)
    plot_fun()
    openxlsx::insertPlot(wb, sheet = sheet_name, width = width, height = height, units = "in", dpi = 300)
  }

  .scale_rows <- function(mat) {
    if (!is.matrix(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
      return(mat)
    }

    scaled <- t(scale(t(mat)))
    scaled[!is.finite(scaled)] <- 0
    scaled
  }

  stats_tbl <- as.data.frame(comparison_result@stats, stringsAsFactors = FALSE)
  if (!"feature_id" %in% colnames(stats_tbl)) {
    stop("`comparison_result$stats` must contain a `feature_id` column.", call. = FALSE)
  }
  if (!"p_value" %in% colnames(stats_tbl) || !"mean_diff" %in% colnames(stats_tbl)) {
    stop("`comparison_result$stats` must contain `p_value` and `mean_diff` columns.", call. = FALSE)
  }
  if (!"mean_group1" %in% colnames(stats_tbl) || !"mean_group2" %in% colnames(stats_tbl)) {
    stop("`comparison_result$stats` must contain `mean_group1` and `mean_group2` columns.", call. = FALSE)
  }

  if (!"feature_label" %in% colnames(stats_tbl)) {
    stats_tbl$feature_label <- stats_tbl$feature_id
  }

  plot_data <- NULL
  if ("plot_data" %in% methods::slotNames(class(comparison_result))) {
    plot_data <- tryCatch(
      comparison_result@plot_data,
      error = function(e) NULL
    )
  }

  response <- NULL
  sample_groups <- NULL
  if (is.list(plot_data) && !is.null(plot_data$response) && !is.null(plot_data$sample_groups)) {
    response <- plot_data$response
    sample_groups <- as.character(plot_data$sample_groups)
    if (!is.matrix(response) || ncol(response) != length(sample_groups)) {
      stop("`comparison_result@plot_data` has inconsistent response matrix dimensions.", call. = FALSE)
    }

    if (is.null(colnames(response))) {
      colnames(response) <- paste0("sample_", seq_len(ncol(response)))
    }
  } else {
    message(
      "`comparison_result` does not contain plot data. Writing the workbook without Boxplot/Scatterplot/Heatmap sheets."
    )
  }

  stats_tbl <- stats_tbl[order(stats_tbl$mean_diff, stats_tbl$p_value, na.last = TRUE), , drop = FALSE]

  group1_display <- .display_label(group1_label)
  group2_display <- .display_label(group2_label)
  comparison_sheet_name <- .make_unique_sheet_name(sprintf("%s_vs_%s", group1_display, group2_display), character())
  readme_title <- sprintf("How to read the %s table:", comparison_sheet_name)

  stats_tbl$feature_label <- as.character(stats_tbl$feature_label)
  missing_feature_label <- is.na(stats_tbl$feature_label) | !nzchar(stats_tbl$feature_label)
  stats_tbl$feature_label[missing_feature_label] <- stats_tbl$feature_id[missing_feature_label]

  group1_mean_col <- paste0("mean_", .clean_label(group1_display))
  group2_mean_col <- paste0("mean_", .clean_label(group2_display))

  export_tbl <- stats_tbl
  export_tbl[[group1_mean_col]] <- export_tbl$mean_group1
  export_tbl[[group2_mean_col]] <- export_tbl$mean_group2
  export_tbl$Difference <- export_tbl$mean_diff
  export_tbl$mean_group1 <- NULL
  export_tbl$mean_group2 <- NULL
  export_tbl$mean_diff <- NULL

  preferred_cols <- c(
    "feature_id",
    "feature_label",
    setdiff(names(export_tbl), c(
      "feature_id",
      "feature_label",
      group1_mean_col,
      group2_mean_col,
      "Difference",
      "n_group1",
      "n_group2",
      "p_value",
      "p_adj",
      "neg_log10_p",
      "significant"
    )),
    group1_mean_col,
    group2_mean_col,
    "Difference",
    "n_group1",
    "n_group2",
    "p_value",
    "p_adj",
    "neg_log10_p",
    "significant"
  )
  preferred_cols <- unique(intersect(preferred_cols, names(export_tbl)))
  export_tbl <- export_tbl[, preferred_cols, drop = FALSE]

  wb <- openxlsx::createWorkbook()
  sheet_names <- character()

  add_sheet_name <- function(base_name) {
    sheet_name <- .make_unique_sheet_name(base_name, sheet_names)
    sheet_names <<- c(sheet_names, sheet_name)
    sheet_name
  }

  readme_sheet <- add_sheet_name("README")
  comparison_sheet <- add_sheet_name(comparison_sheet_name)

  openxlsx::addWorksheet(wb, readme_sheet, gridLines = FALSE)
  openxlsx::addWorksheet(wb, comparison_sheet, gridLines = FALSE)
  readme <- c(
    readme_title,
    "",
    sprintf("PRISM observed: %s (n=%d) vs. %s (n=%d)", group1_display, comparison_result@group1_n_models, group2_display, comparison_result@group2_n_models),
    "",
    "Each row is one PRISM compound.",
    sprintf("Drugs are sorted by Difference = mean %s response - mean %s response.", group1_display, group2_display),
    sprintf("The most negative values at the top indicate drugs to which %s cell lines are most sensitive relative to %s.", group1_display, group2_display),
    "Lower PRISM response values indicate greater growth inhibition, meaning greater sensitivity.",
    "",
    "Key columns:",
    "feature_id   = treatment or drug identifier",
    "feature_label = label used for plotting when available",
    sprintf("%s = mean %s response", group1_mean_col, group1_display),
    sprintf("%s = mean %s response", group2_mean_col, group2_display),
    "Difference   = mean(group1) - mean(group2)",
    "p_value      = Welch t-test p-value",
    "p_adj        = adjusted p-value",
    "significant  = p_adj < 0.05 and |mean_diff| above the threshold"
  )
  openxlsx::writeData(wb, readme_sheet, readme, startCol = 1, startRow = 1)
  openxlsx::setColWidths(wb, readme_sheet, cols = 1, widths = 120.71)

  openxlsx::writeData(wb, comparison_sheet, export_tbl, startRow = 1, rowNames = FALSE)
  openxlsx::freezePane(wb, comparison_sheet, firstActiveRow = 2)
  openxlsx::addFilter(wb, comparison_sheet, row = 1, cols = seq_len(ncol(export_tbl)))
  openxlsx::setColWidths(wb, comparison_sheet, cols = seq_len(ncol(export_tbl)), widths = "auto")

  plot_stats <- stats_tbl
  plot_stats$neg_log10_p <- -log10(plot_stats$p_value)
  point_colors <- ifelse(plot_stats$significant, "#C0392B", "gray70")
  volcano_labels <- plot_stats[order(plot_stats$p_value, -abs(plot_stats$mean_diff)), , drop = FALSE]
  volcano_labels <- utils::head(volcano_labels, n_label_volcano)

  if (include_volcano) {
    .make_plot_sheet(wb, add_sheet_name("Volcano"), function() {
      graphics::par(mar = c(5, 5, 4, 2) + 0.1)
      graphics::plot(
        plot_stats$mean_diff,
        plot_stats$neg_log10_p,
        pch = 16,
        col = point_colors,
        xlab = sprintf("Mean difference (%s - %s)", group1_display, group2_display),
        ylab = expression(-log[10](p-value)),
        main = sprintf("Drug sensitivity: %s vs %s", group1_display, group2_display)
      )
      graphics::abline(h = -log10(0.05), lty = 2, col = "gray55")
      graphics::abline(v = 0, lty = 2, col = "gray55")
      if (nrow(volcano_labels) > 0L) {
        graphics::text(
          volcano_labels$mean_diff,
          volcano_labels$neg_log10_p,
          labels = volcano_labels$feature_label,
          pos = 3,
          cex = 0.7,
          xpd = NA
        )
      }
    })
  } else {
    openxlsx::addWorksheet(wb, add_sheet_name("Volcano"), gridLines = FALSE)
  }

  if (!is.null(response) && !is.null(sample_groups)) {
    top_negative <- plot_stats[which.min(plot_stats$mean_diff), , drop = FALSE]
    top_positive <- plot_stats[which.max(plot_stats$mean_diff), , drop = FALSE]

    .make_plot_sheet(wb, add_sheet_name("Boxplot1"), function() {
      feature_row <- top_negative[1, , drop = FALSE]
      feature_id <- feature_row$feature_id
      feature_values <- as.numeric(response[feature_id, ])
      valid_idx <- is.finite(feature_values) & !is.na(sample_groups)
      plot_groups <- factor(sample_groups[valid_idx], levels = c(group1_display, group2_display))
      plot_values <- feature_values[valid_idx]

      graphics::par(mar = c(6, 5, 4, 2) + 0.1)
      graphics::boxplot(
        plot_values ~ plot_groups,
        col = c("#4E79A7", "#E15759"),
        outline = FALSE,
        ylab = sprintf("Response (%s)", comparison_result@unit),
        main = sprintf("Most sensitive in %s: %s", group1_display, .feature_title(feature_row))
      )
      graphics::stripchart(
        plot_values ~ plot_groups,
        vertical = TRUE,
        method = "jitter",
        add = TRUE,
        pch = 16,
        col = grDevices::adjustcolor("gray20", alpha.f = 0.45)
      )
    }, width = 8.5, height = 5.5)

    .make_plot_sheet(wb, add_sheet_name("Boxplot2"), function() {
      feature_row <- top_positive[1, , drop = FALSE]
      feature_id <- feature_row$feature_id
      feature_values <- as.numeric(response[feature_id, ])
      valid_idx <- is.finite(feature_values) & !is.na(sample_groups)
      plot_groups <- factor(sample_groups[valid_idx], levels = c(group1_display, group2_display))
      plot_values <- feature_values[valid_idx]

      graphics::par(mar = c(6, 5, 4, 2) + 0.1)
      graphics::boxplot(
        plot_values ~ plot_groups,
        col = c("#4E79A7", "#E15759"),
        outline = FALSE,
        ylab = sprintf("Response (%s)", comparison_result@unit),
        main = sprintf("Most sensitive in %s: %s", group2_display, .feature_title(feature_row))
      )
      graphics::stripchart(
        plot_values ~ plot_groups,
        vertical = TRUE,
        method = "jitter",
        add = TRUE,
        pch = 16,
        col = grDevices::adjustcolor("gray20", alpha.f = 0.45)
      )
    }, width = 8.5, height = 5.5)

    .make_plot_sheet(wb, add_sheet_name("Scatterplot"), function() {
      graphics::par(mar = c(5, 5, 4, 2) + 0.1)
      graphics::plot(
        plot_stats$mean_group1,
        plot_stats$mean_group2,
        pch = 16,
        col = point_colors,
        xlab = sprintf("Mean response in %s", group1_display),
        ylab = sprintf("Mean response in %s", group2_display),
        main = sprintf("Mean response: %s vs %s", group1_display, group2_display)
      )
      graphics::abline(0, 1, lty = 2, col = "gray55")
      if (nrow(volcano_labels) > 0L) {
        graphics::text(
          volcano_labels$mean_group1,
          volcano_labels$mean_group2,
          labels = volcano_labels$feature_label,
          pos = 3,
          cex = 0.7,
          xpd = NA
        )
      }
    })

    .make_plot_sheet(wb, add_sheet_name("Heatmap"), function() {
      heat_feature_count <- min(25L, nrow(response))
      heat_features <- plot_stats$feature_id[order(plot_stats$p_value, -abs(plot_stats$mean_diff))][seq_len(heat_feature_count)]
      heat_labels <- plot_stats$feature_label[match(heat_features, plot_stats$feature_id)]
      heat_values <- response[heat_features, , drop = FALSE]

      group_order <- order(factor(sample_groups, levels = c(group1_display, group2_display)), colnames(response))
      heat_values <- heat_values[, group_order, drop = FALSE]
      heat_matrix <- .scale_rows(as.matrix(heat_values))
      heat_matrix <- heat_matrix[rev(seq_len(nrow(heat_matrix))), , drop = FALSE]

      graphics::par(mar = c(8, 9, 4, 2) + 0.1)
      graphics::image(
        x = seq_len(ncol(heat_matrix)),
        y = seq_len(nrow(heat_matrix)),
        z = t(heat_matrix),
        col = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        axes = FALSE,
        xlab = "Samples",
        ylab = "Top features by significance",
        main = sprintf("Heatmap: %s vs %s", group1_display, group2_display)
      )
      graphics::axis(1, at = seq_len(ncol(heat_matrix)), labels = colnames(heat_values), las = 2, cex.axis = 0.5)
      graphics::axis(2, at = seq_len(nrow(heat_matrix)), labels = rev(heat_labels), las = 1, cex.axis = 0.6)
      graphics::box()
    }, width = 8.5, height = 6.5)
  } else {
    openxlsx::addWorksheet(wb, add_sheet_name("Boxplot1"), gridLines = FALSE)
    openxlsx::addWorksheet(wb, add_sheet_name("Boxplot2"), gridLines = FALSE)
    openxlsx::addWorksheet(wb, add_sheet_name("Scatterplot"), gridLines = FALSE)
    openxlsx::addWorksheet(wb, add_sheet_name("Heatmap"), gridLines = FALSE)
  }

  openxlsx::saveWorkbook(wb, file = output_file, overwrite = overwrite)
  invisible(output_file)
}

#' Export comparison results
#'
#' @param object A `DrugSensitivityComparison` object.
#' @param ... Additional arguments passed to the export method.
#' @keywords internal
methods::setGeneric(
  name = "export",
  def = function(object, ...) {
    standardGeneric("export")
  }
)
