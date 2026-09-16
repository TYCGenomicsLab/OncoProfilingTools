library(readr)
library(dplyr)
library(ggplot2)
library(enrichR)

setEnrichrSite("Enrichr")

source("R/agent-chea.R")

input_file <- "output/cms4_fc2/cms4_fc2_degs.csv"
output_dir <- "output/chea_cms4"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Reading CMS4 DEG list...\n")

deg <- read_csv(
  input_file,
  show_col_types = FALSE
)

possible_gene_columns <- c(
  "gene_symbol",
  "Gene",
  "gene",
  "SYMBOL",
  "Geneid"
)

gene_column <- possible_gene_columns[
  possible_gene_columns %in% names(deg)
][1]

if (is.na(gene_column)) {
  stop(
    paste(
      "No gene-symbol column found. Columns:",
      paste(names(deg), collapse = ", ")
    ),
    call. = FALSE
  )
}

genes <- unique(na.omit(deg[[gene_column]]))
genes <- genes[genes != ""]

cat("Using gene column:", gene_column, "\n")
cat("Input genes:", length(genes), "\n")
cat("Running ChEA enrichment...\n")

result <- run_chea_agent(
  genes = genes,
  database = "ChEA_2022"
)

cat(result$summary, "\n")

write_csv(
  result$results,
  file.path(output_dir, "chea_results.csv")
)

if (nrow(result$results) > 0) {
  plot_df <- result$results %>%
    filter(
      !is.na(Adjusted.P.value),
      !is.na(Combined.Score)
    ) %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = 20) %>%
    mutate(
      Term = factor(Term, levels = rev(Term)),
      neg_log10_padj = -log10(
        pmax(Adjusted.P.value, .Machine$double.xmin)
      )
    )

  p <- ggplot(
    plot_df,
    aes(
      x = Combined.Score,
      y = Term,
      size = neg_log10_padj,
      color = Adjusted.P.value
    )
  ) +
    geom_point(alpha = 0.9) +
    labs(
      title = "Top ChEA Transcription Factor Enrichment",
      subtitle = "CMS4 genes analyzed using ChEA 2022",
      x = "Combined Score",
      y = NULL,
      size = "-log10 adjusted p-value",
      color = "Adjusted p-value"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    file.path(output_dir, "chea_tf_dotplot.png"),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )
}

cat("ChEA analysis completed successfully.\n")