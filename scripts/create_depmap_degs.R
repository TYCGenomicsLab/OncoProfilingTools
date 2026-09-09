model <- read.csv("data/Model.csv")
expr <- read.csv("data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                 check.names = FALSE)

bowel_ids <- model$ModelID[model$OncotreeLineage == "Bowel"]

expr_bowel <- expr[expr[[1]] %in% bowel_ids, ]
expr_other <- expr[!(expr[[1]] %in% bowel_ids), ]

genes_raw <- names(expr)[-1]

clean_gene <- function(x) {
  x <- gsub(" \\([0-9]+\\)$", "", x)
  x <- gsub("\\.\\.\\.[0-9]+$", "", x)
  x
}

deg <- data.frame(
  gene_symbol = clean_gene(genes_raw),
  mean_bowel = sapply(genes_raw, function(g) mean(expr_bowel[[g]], na.rm = TRUE)),
  mean_other = sapply(genes_raw, function(g) mean(expr_other[[g]], na.rm = TRUE))
)

deg$logFC <- deg$mean_bowel - deg$mean_other
deg <- deg[order(abs(deg$logFC), decreasing = TRUE), ]

dir.create("output/depmap_bowel", recursive = TRUE, showWarnings = FALSE)

write.csv(head(deg, 100),
          "output/depmap_bowel/depmap_bowel_top100_degs.csv",
          row.names = FALSE)

message("Done: output/depmap_bowel/depmap_bowel_top100_degs.csv")
