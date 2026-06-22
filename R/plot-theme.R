# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

#' OncoProfilingTools ggplot theme
#'
#' A minimal theme for publication-style oncology profiling figures.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A ggplot2 theme.
#'
#' @export
theme_onco <- function(
  base_size = 12,
  base_family = "sans"
) {
  ggplot2::theme_minimal(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(),
      legend.title = ggplot2::element_text(face = "bold")
    )
}
