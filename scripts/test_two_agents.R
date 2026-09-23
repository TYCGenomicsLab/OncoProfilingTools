library(readr)

source("R/agent-report.R")
source("R/agent-kegg.R")
source("R/agent-go.R")

dir.create("output", showWarnings = FALSE)

deg_data <- readr::read_csv("inst/extdata/sample_deg.csv", show_col_types = FALSE)
genes <- unique(deg_data$gene_symbol)

kegg_result <- run_kegg_agent(genes)
go_result <- run_go_agent(genes)

write.csv(kegg_result$results, "output/kegg_results.csv", row.names = FALSE)
write.csv(go_result$results, "output/go_results.csv", row.names = FALSE)

generate_two_agent_report(
  kegg_result,
  go_result,
  output_file = "output/combined_report.md"
)

message("Done. Check output folder.")
