source("R/agent-gsva.R")

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

set.seed(123)

genes <- paste0("GENE", 1:100)

expression_matrix <- matrix(
  rnorm(100 * 6),
  nrow = 100,
  ncol = 6,
  dimnames = list(
    genes,
    paste0("Sample", 1:6)
  )
)

gene_sets <- list(
  PATHWAY_A = genes[1:20],
  PATHWAY_B = genes[21:45],
  PATHWAY_C = genes[46:75]
)

result <- run_gsva_agent(
  expression_matrix = expression_matrix,
  gene_sets = gene_sets,
  kcdf = "Gaussian"
)

dir.create(
  "output/gsva_test",
  recursive = TRUE,
  showWarnings = FALSE
)

png(
  "output/gsva_test/gsva_heatmap.png",
  width = 1600,
  height = 1000,
  res = 180
)

pheatmap::pheatmap(
  result$results,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "GSVA Pathway Activity Scores",
  border_color = NA
)

dev.off()

cat("GSVA heatmap generated successfully.\n")