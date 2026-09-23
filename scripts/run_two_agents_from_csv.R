source("R/agent-report.R")
source("R/agent-kegg.R")
source("R/agent-go.R")

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) args[1] else "inst/extdata/sample_deg.csv"
output_dir <- if (length(args) >= 2) args[2] else "output"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

deg_data <- readr::read_csv(input_file, show_col_types = FALSE)

if (!"gene_symbol" %in% names(deg_data)) {
  stop("Input CSV must contain a column named 'gene_symbol'.")
}

genes <- unique(deg_data$gene_symbol)

kegg_result <- run_kegg_agent(genes)
go_result <- run_go_agent(genes)

write.csv(kegg_result$results, file.path(output_dir, "kegg_results.csv"), row.names = FALSE)
write.csv(go_result$results, file.path(output_dir, "go_results.csv"), row.names = FALSE)

generate_two_agent_report(
  kegg_result,
  go_result,
  output_file = file.path(output_dir, "combined_report.md")
)

message("Two-agent pipeline completed successfully.")
