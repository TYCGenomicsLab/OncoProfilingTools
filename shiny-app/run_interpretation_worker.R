args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Expected input RDS, output RDS, and interpretation helper paths.", call. = FALSE)
}

input_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_path <- args[[2L]]
helper_path <- normalizePath(args[[3L]], mustWork = TRUE)

source(helper_path, local = TRUE)
payload <- readRDS(input_path)
bundle <- generate_provider_interpretation_bundle(payload$data_by_agent, payload$settings)

temporary_output <- paste0(output_path, ".tmp")
saveRDS(bundle, temporary_output)
if (!file.rename(temporary_output, output_path)) {
  stop("Could not publish the local interpretation result.", call. = FALSE)
}
