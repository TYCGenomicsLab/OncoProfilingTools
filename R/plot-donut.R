# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Create a donut plot
#'
#' Creates a donut plot from assay-level metadata stored in an `OncoExperiment`.
#'
#' @param object An `OncoExperiment` object.
#' @param variable Metadata field to plot. Can be unquoted, a character string,
#'   or an enum value such as `ModelMetadataFields$OncotreeLineage`.
#' @param assay A single character string specifying which assay experiment to
#'   use. For example, `"Expression"` or `"PRISM"`.
#' @param variable_expr Captured unevaluated variable expression. Internal.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_donut(exp, ModelMetadataFields$OncotreeLineage, assay = "Expression")
#' plot_donut(exp, "Sex", assay = "Expression")
#' plot_donut(exp, OncotreeLineage, assay = "Expression")
#' }
#'
#' @export
methods::setGeneric(
  name = "plot_donut",
  def = function(
    object,
    variable,
    assay,
    variable_expr = substitute(variable),
    ...
  ) {
    standardGeneric("plot_donut")
  },
  signature = "object"
)

#' Create a donut plot from model metadata
#'
#' Generates a donut plot that shows the distribution of categorical metadata
#' fields from a selected assay in an `OncoExperiment`.
#'
#' @param object An `OncoExperiment` object.
#' @param variable Metadata field to plot.
#' @param assay Assay name to use.
#' @param variable_expr Captured unevaluated variable expression. Internal.
#' @param show_na Logical. Whether missing values should be shown as `"Missing"`.
#' @param ... Additional arguments passed to `.apply_onco_plot_theme()`.
#'
#' @return A `ggplot` object.
#'
#' @include plot-theme.R
#' @include enum-ModelMetadataFields.R
#' @export
methods::setMethod(
  f = "plot_donut",
  signature = signature(object = "OncoExperiment"),
  definition = function(
    object,
    variable,
    assay,
    variable_expr = substitute(variable),
    show_na = FALSE,
    ...
  ) {
    ## ---------------------------------------------------
    ## Validate arguments
    ## ---------------------------------------------------

    if (missing(variable)) {
      stop(
        "Please specify a metadata field to plot. Call `ModelMetadataFields` ",
        "for available ExpressionAssay metadata fields.",
        call. = FALSE
      )
    }

    if (missing(assay)) {
      stop("Please specify an assay experiment to use.", call. = FALSE)
    }

    assay_name <- .resolve_plot_donut_assay(assay)

    if (!is.logical(show_na) || length(show_na) != 1L || is.na(show_na)) {
      stop("`show_na` must be TRUE or FALSE.", call. = FALSE)
    }

    ## ---------------------------------------------------
    ## Extract assay experiment
    ## ---------------------------------------------------

    experiments <- MultiAssayExperiment::experiments(object)
    experiment_names <- names(experiments)

    if (is.null(experiment_names) || length(experiment_names) == 0L) {
      stop(
        "This OncoExperiment does not contain any assay experiments.",
        call. = FALSE
      )
    }

    if (!assay_name %in% experiment_names) {
      stop(
        "Assay experiment `",
        assay_name,
        "` not found in this OncoExperiment.\n",
        "Available assays: ",
        paste(experiment_names, collapse = ", "),
        call. = FALSE
      )
    }

    assay_experiment <- experiments[[assay_name]]

    if (!inherits(assay_experiment, "SummarizedExperiment")) {
      stop(
        "Assay experiment `",
        assay_name,
        "` is not a SummarizedExperiment. ",
        "Cannot extract `colData()` for plotting.",
        call. = FALSE
      )
    }

    assay_coldata <- as.data.frame(
      SummarizedExperiment::colData(assay_experiment),
      stringsAsFactors = FALSE
    )

    if (ncol(assay_coldata) == 0L) {
      stop(
        "Assay experiment `",
        assay_name,
        "` does not contain any column metadata.",
        call. = FALSE
      )
    }

    ## ---------------------------------------------------
    ## Resolve variable
    ## ---------------------------------------------------

    variable_name <- .resolve_plot_donut_variable(
      variable = variable,
      variable_expr = variable_expr
    )

    if (!variable_name %in% colnames(assay_coldata)) {
      stop(
        "Variable `",
        variable_name,
        "` not found in column metadata of assay `",
        assay_name,
        "`.\nAvailable variables: ",
        paste(colnames(assay_coldata), collapse = ", "),
        call. = FALSE
      )
    }

    ## ---------------------------------------------------
    ## Build summary data frame for plotting
    ## ---------------------------------------------------

    values <- assay_coldata[[variable_name]]

    plot_df <- data.frame(
      variable = as.character(values),
      stringsAsFactors = FALSE
    )

    is_missing <- is.na(plot_df$variable) | !nzchar(plot_df$variable)

    if (isTRUE(show_na)) {
      plot_df$variable[is_missing] <- "Missing"
    } else {
      plot_df <- plot_df[!is_missing, , drop = FALSE]
    }

    if (nrow(plot_df) == 0L) {
      stop(
        "No non-missing values available for variable `",
        variable_name,
        "` in assay `",
        assay_name,
        "`.",
        call. = FALSE
      )
    }

    plot_df <- plot_df |>
      dplyr::count(.data$variable, name = "count") |>
      dplyr::arrange(dplyr::desc(.data$count), .data$variable) |>
      dplyr::mutate(
        fraction = .data$count / sum(.data$count),
        percent = .data$fraction * 100,
        ymax = cumsum(.data$fraction),
        ymin = dplyr::lag(.data$ymax, default = 0),
        label_position = (.data$ymax + .data$ymin) / 2,
        label = as.character(.data$count)
      )

    ## ---------------------------------------------------
    ## Create donut plot
    ## ---------------------------------------------------

    ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        ymax = .data$ymax,
        ymin = .data$ymin,
        xmax = 4,
        xmin = 2.4,
        fill = .data$variable
      )
    ) +
      ggplot2::geom_rect(
        color = "white",
        linewidth = 0.4
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          x = 3.2,
          y = .data$label_position,
          label = .data$label
        ),
        inherit.aes = FALSE,
        size = 3
      ) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::xlim(0, 4) +
      ggplot2::labs(
        title = variable_name,
        subtitle = paste0("n = ", sum(plot_df$count)),
        fill = NULL
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5
        ),
        legend.position = "right",
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }
)

#' Resolve assay argument for plot_donut
#'
#' @param assay Assay name or enum value.
#'
#' @return A single character string.
#'
#' @keywords internal
.resolve_plot_donut_assay <- function(assay) {
  assay_name <- as.character(assay)

  if (
    !is.character(assay_name) ||
      length(assay_name) != 1L ||
      is.na(assay_name) ||
      !nzchar(assay_name)
  ) {
    stop(
      "`assay` must resolve to a single non-empty character string.",
      call. = FALSE
    )
  }

  assay_name
}

#' Resolve variable argument for plot_donut
#'
#' @param variable Evaluated variable argument.
#' @param variable_expr Captured unevaluated variable expression.
#'
#' @return A single character string.
#'
#' @keywords internal
.resolve_plot_donut_variable <- function(variable, variable_expr) {
  ## ---------------------------------------------------
  ## Case 1: character string or enum value
  ## ---------------------------------------------------

  variable_value <- tryCatch(
    variable,
    error = function(e) NULL
  )

  if (!is.null(variable_value)) {
    variable_name <- as.character(variable_value)

    if (
      is.character(variable_name) &&
        length(variable_name) == 1L &&
        !is.na(variable_name) &&
        nzchar(variable_name)
    ) {
      return(variable_name)
    }
  }

  ## ---------------------------------------------------
  ## Case 2: unquoted column name
  ## ---------------------------------------------------

  if (is.symbol(variable_expr)) {
    variable_name <- as.character(variable_expr)

    if (
      is.character(variable_name) &&
        length(variable_name) == 1L &&
        !is.na(variable_name) &&
        nzchar(variable_name) &&
        variable_name != "variable"
    ) {
      return(variable_name)
    }
  }

  stop(
    "`variable` must be an unquoted column name, a character string, ",
    "or an enum value.",
    call. = FALSE
  )
}

#' Acquire the plotting theme from the environment
#'
#' If `theme_onco()` is not available, returns a generic ggplot minimal theme.
#'
#' @param base_size Base font size.
#'
#' @return A ggplot theme.
#'
#' @keywords internal
.apply_onco_plot_theme <- function(base_size = 12) {
  if (exists("theme_onco", mode = "function", inherits = TRUE)) {
    return(
      get("theme_onco", mode = "function", inherits = TRUE)(
        base_size = base_size
      )
    )
  }

  ggplot2::theme_minimal(base_size = base_size)
}
