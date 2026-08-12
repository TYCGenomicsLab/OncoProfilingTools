library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)

deg_data <- readr::read_csv("inst/extdata/sample_deg.csv", show_col_types = FALSE)

genes <- unique(deg_data$gene_symbol)

mapped_genes <- clusterProfiler::bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)

print(mapped_genes)
