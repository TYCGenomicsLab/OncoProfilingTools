source("R/agent-gsva.R")

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

print(result$summary)
print(dim(result$results))
print(result$results)

dir.create(
  "output/gsva_test",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  result$results,
  "output/gsva_test/gsva_scores.csv",
  row.names = TRUE
)

cat("GSVA test completed successfully.\n")