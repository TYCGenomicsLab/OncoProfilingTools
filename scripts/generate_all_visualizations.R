#!/usr/bin/env Rscript

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  stringsAsFactors = FALSE
)

required_packages <- c(
  "readr",
  "dplyr",
  "ggplot2",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------

output_directory <- "output/visualizations"

go_input <- "output/cms4/go_results.csv"
kegg_input <- "output/cms4/kegg_results.csv"

go_plot_path <- file.path(
  output_directory,
  "go_biological_process_dotplot.png"
)

kegg_plot_path <- file.path(
  output_directory,
  "kegg_pathway_dotplot.png"
)

gsva_plot_path <- "output/gsva_bowel/gsva_hallmark_heatmap.png"
chea_plot_path <- "output/chea_cms4/chea_tf_dotplot.png"

manifest_path <- file.path(
  output_directory,
  "visualization_manifest.csv"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

remove_stale_file <- function(path) {
  if (file.exists(path)) {
    unlink(path)
    message("Removed stale visualization: ", path)
  }
}

prepare_enrichment_data <- function(
  data,
  maximum_terms = 15
) {
  required_columns <- c(
    "Description",
    "p.adjust",
    "Count"
  )

  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  cleaned <- data %>%
    mutate(
      Description = trimws(as.character(Description)),
      p.adjust = safe_numeric(p.adjust),
      Count = safe_numeric(Count),
      enrichment_score = ifelse(
        !is.na(p.adjust) & p.adjust > 0,
        -log10(p.adjust),
        NA_real_
      )
    ) %>%
    filter(
      !is.na(Description),
      Description != "",
      !is.na(p.adjust),
      p.adjust > 0,
      p.adjust <= 1,
      is.finite(enrichment_score),
      !is.na(Count),
      Count > 0
    ) %>%
    arrange(
      p.adjust,
      desc(Count)
    ) %>%
    distinct(
      Description,
      .keep_all = TRUE
    ) %>%
    slice_head(n = maximum_terms)

  cleaned
}

create_dotplot <- function(
  data,
  title,
  subtitle,
  output_path
) {
  if (nrow(data) == 0) {
    stop("Cannot create a dot plot from zero valid rows.")
  }

  plot_data <- data %>%
    mutate(
      Description = factor(
        Description,
        levels = rev(Description)
      )
    )

  plot_object <- ggplot(
    plot_data,
    aes(
      x = enrichment_score,
      y = Description,
      size = Count,
      color = enrichment_score
    )
  ) +
    geom_point(
      alpha = 0.9
    ) +
    scale_size_continuous(
      name = "Gene count",
      range = c(4, 13),
      breaks = pretty_breaks(n = 5)
    ) +
    scale_color_viridis_c(
      name = expression(-log[10]("adjusted p-value")),
      option = "plasma",
      direction = 1
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.03, 0.12))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      caption = paste(
        "Larger values indicate stronger enrichment.",
        "Only valid terms with adjusted p-values greater than zero are shown."
      )
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 22,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        size = 14,
        margin = margin(b = 20)
      ),
      plot.caption = element_text(
        color = "grey40",
        hjust = 0,
        margin = margin(t = 16)
      ),
      axis.title.x = element_text(
        size = 14,
        margin = margin(t = 12)
      ),
      axis.text.x = element_text(
        size = 11
      ),
      axis.text.y = element_text(
        size = 11
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      plot.margin = margin(
        t = 25,
        r = 30,
        b = 25,
        l = 25
      )
    )

  ggsave(
    filename = output_path,
    plot = plot_object,
    width = 14,
    height = 9,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  invisible(output_path)
}

# -------------------------------------------------------------------
# GO Biological Process visualization
# -------------------------------------------------------------------

cat("\nReading GO results...\n")

go_ready <- FALSE
go_message <- "Not processed"

if (!file.exists(go_input)) {
  go_message <- paste("Input file does not exist:", go_input)
  remove_stale_file(go_plot_path)
  message(go_message)
} else {
  go_results <- read_csv(
    go_input,
    show_col_types = FALSE
  )

  cat("GO rows:", nrow(go_results), "\n")
  cat(
    "GO columns:",
    paste(names(go_results), collapse = ", "),
    "\n"
  )

  if (nrow(go_results) == 0) {
    go_message <- "GO result file contains zero rows"
    remove_stale_file(go_plot_path)
    message(go_message)
  } else {
    go_clean <- prepare_enrichment_data(
      go_results,
      maximum_terms = 15
    )

    cat("Valid GO rows for plotting:", nrow(go_clean), "\n")

    if (nrow(go_clean) == 0) {
      go_message <- "GO result file contains no valid plottable rows"
      remove_stale_file(go_plot_path)
      message(go_message)
    } else {
      create_dotplot(
        data = go_clean,
        title = "GO Biological Process Enrichment",
        subtitle = paste(
          "Top biological processes identified",
          "from the CMS4 gene set"
        ),
        output_path = go_plot_path
      )

      go_ready <- file.exists(go_plot_path)
      go_message <- paste(
        nrow(go_clean),
        "GO terms plotted successfully"
      )

      message("GO visualization generated successfully.")
    }
  }
}

# -------------------------------------------------------------------
# KEGG visualization
# -------------------------------------------------------------------

cat("\nReading KEGG results...\n")

kegg_ready <- FALSE
kegg_message <- "Not processed"

if (!file.exists(kegg_input)) {
  kegg_message <- paste("Input file does not exist:", kegg_input)
  remove_stale_file(kegg_plot_path)
  message(kegg_message)
} else {
  kegg_results <- read_csv(
    kegg_input,
    show_col_types = FALSE
  )

  cat("KEGG rows:", nrow(kegg_results), "\n")
  cat(
    "KEGG columns:",
    paste(names(kegg_results), collapse = ", "),
    "\n"
  )

  if (nrow(kegg_results) == 0) {
    kegg_message <- paste(
      "KEGG analysis returned zero enriched pathways;",
      "no KEGG plot was generated"
    )

    remove_stale_file(kegg_plot_path)
    message(kegg_message)
  } else {
    kegg_clean <- prepare_enrichment_data(
      kegg_results,
      maximum_terms = 15
    )

    cat(
      "Valid KEGG rows for plotting:",
      nrow(kegg_clean),
      "\n"
    )

    if (nrow(kegg_clean) == 0) {
      kegg_message <- paste(
        "KEGG results contained no valid plottable pathways;",
        "no KEGG plot was generated"
      )

      remove_stale_file(kegg_plot_path)
      message(kegg_message)
    } else {
      create_dotplot(
        data = kegg_clean,
        title = "KEGG Pathway Enrichment",
        subtitle = paste(
          "Top signaling and disease pathways identified",
          "from the CMS4 gene set"
        ),
        output_path = kegg_plot_path
      )

      kegg_ready <- file.exists(kegg_plot_path)
      kegg_message <- paste(
        nrow(kegg_clean),
        "KEGG pathways plotted successfully"
      )

      message("KEGG visualization generated successfully.")
    }
  }
}

# -------------------------------------------------------------------
# Existing GSVA and ChEA visualizations
# -------------------------------------------------------------------

gsva_ready <- file.exists(gsva_plot_path)
chea_ready <- file.exists(chea_plot_path)

gsva_message <- if (gsva_ready) {
  "GSVA Hallmark heatmap exists"
} else {
  paste("GSVA visualization is missing:", gsva_plot_path)
}

chea_message <- if (chea_ready) {
  "ChEA transcription-factor dot plot exists"
} else {
  paste("ChEA visualization is missing:", chea_plot_path)
}

# -------------------------------------------------------------------
# Manifest
# -------------------------------------------------------------------

manifest <- tibble::tibble(
  agent = c(
    "GO",
    "KEGG",
    "GSVA",
    "ChEA"
  ),
  visualization = c(
    go_plot_path,
    kegg_plot_path,
    gsva_plot_path,
    chea_plot_path
  ),
  ready = c(
    go_ready,
    kegg_ready,
    gsva_ready,
    chea_ready
  ),
  status = c(
    go_message,
    kegg_message,
    gsva_message,
    chea_message
  )
)

write_csv(
  manifest,
  manifest_path
)

# -------------------------------------------------------------------
# Final status
# -------------------------------------------------------------------

cat("\nVisualization status\n")
cat("--------------------\n")
cat("GO:   ", go_ready, "\n")
cat("KEGG: ", kegg_ready, "\n")
cat("GSVA: ", gsva_ready, "\n")
cat("ChEA: ", chea_ready, "\n")

cat("\nDetails\n")
cat("-------\n")

for (index in seq_len(nrow(manifest))) {
  cat(
    manifest$agent[index],
    ": ",
    manifest$status[index],
    "\n",
    sep = ""
  )
}

cat(
  "\nManifest saved to:",
  manifest_path,
  "\n"
)

if (all(manifest$ready)) {
  cat("\nAll four agent visualizations are ready.\n")
} else {
  cat(
    "\n",
    sum(manifest$ready),
    "of",
    nrow(manifest),
    "agent visualizations are currently ready.\n"
  )
}

if (!kegg_ready && file.exists(kegg_input)) {
  cat(
    paste0(
      "\nKEGG is correctly marked FALSE because ",
      "the current KEGG results file contains no ",
      "valid enriched pathways.\n"
    )
  )
}
