library(readr)
library(dplyr)
library(msigdbr)
library(pheatmap)

source("R/agent-gsva.R")

input_file <- "data/OmicsExpression_Bowel_TPMLogp1.csv"
output_dir <- "output/gsva_bowel"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading bowel expression matrix...\n")

expr_df <- read_csv(
  input_file,
  show_col_types = FALSE
)

sample_ids <- expr_df[[1]]

expr_values <- as.matrix(expr_df[, -1])

storage.mode(expr_values) <- "numeric"

gene_names_raw <- colnames(expr_df)[-1]

gene_symbols <- sub(
  "\\s*\\([^)]*\\)$",
  "",
  gene_names_raw
)

colnames(expr_values) <- gene_symbols
rownames(expr_values) <- sample_ids

cat(
  "Loaded",
  nrow(expr_values),
  "models and",
  ncol(expr_values),
  "genes.\n"
)

cat("Transposing to genes x models...\n")

expression_matrix <- t(expr_values)

expression_matrix <- expression_matrix[
  !duplicated(rownames(expression_matrix)),
  ,
  drop = FALSE
]

cat("Loading Hallmark gene sets...\n")

hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_sets <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

cat("Running GSVA...\n")

result <- run_gsva_agent(
  expression_matrix = expression_matrix,
  gene_sets = hallmark_sets,
  kcdf = "Gaussian",
  min_size = 10,
  max_size = 500
)

write.csv(
  result$results,
  file.path(output_dir, "gsva_hallmark_scores.csv"),
  row.names = TRUE
)

cat(result$summary, "\n")

cat("Generating heatmap...\n")

top_n <- min(20, nrow(result$results))

pathway_variance <- apply(
  result$results,
  1,
  var,
  na.rm = TRUE
)

top_pathways <- names(
  sort(pathway_variance, decreasing = TRUE)
)[seq_len(top_n)]

heatmap_matrix <- result$results[
  top_pathways,
  ,
  drop = FALSE
]

png(
  file.path(output_dir, "gsva_hallmark_heatmap.png"),
  width = 2400,
  height = 1800,
  res = 220
)

pheatmap(
  heatmap_matrix,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  main = "Top Variable Hallmark Pathway Activity in Bowel Models",
  border_color = NA
)

dev.off()

cat("Real GSVA bowel analysis completed successfully.\n")