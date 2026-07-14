suppressPackageStartupMessages({
  library(readr)
})

`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0) {
    fallback
  } else {
    value
  }
}

project_root <- normalizePath(".", mustWork = TRUE)

worker_script <- file.path(
  project_root,
  "shiny-app",
  "run_agent_worker.R"
)

test_directory <- file.path(
  project_root,
  "shiny-app",
  "runtime-test"
)

dir.create(
  test_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

deg_file <- file.path(
  project_root,
  "output",
  "cms4_fc2",
  "cms4_fc2_degs.csv"
)

expression_file <- file.path(
  project_root,
  "data",
  "OmicsExpression_Bowel_TPMLogp1.csv"
)

if (!file.exists(worker_script)) {
  stop("Worker script is missing: ", worker_script)
}

if (!file.exists(deg_file)) {
  stop("DEG input file is missing: ", deg_file)
}

if (!file.exists(expression_file)) {
  stop("Expression input file is missing: ", expression_file)
}

cat("=== PREPARING REAL TEST INPUTS ===\n")

deg_data <- read_csv(
  deg_file,
  show_col_types = FALSE,
  name_repair = "unique"
)

expression_data <- read_csv(
  expression_file,
  show_col_types = FALSE,
  name_repair = "unique"
)

deg_input_rds <- file.path(
  test_directory,
  "deg-input.rds"
)

expression_input_rds <- file.path(
  test_directory,
  "expression-input.rds"
)

saveRDS(
  deg_data,
  deg_input_rds
)

saveRDS(
  expression_data,
  expression_input_rds
)

cat(
  "DEG dimensions:",
  nrow(deg_data),
  "x",
  ncol(deg_data),
  "\n"
)

cat(
  "Expression dimensions:",
  nrow(expression_data),
  "x",
  ncol(expression_data),
  "\n"
)

agent_inputs <- list(
  go = deg_input_rds,
  kegg = deg_input_rds,
  chea = deg_input_rds,
  gsva = expression_input_rds
)

validation <- list()

for (agent_name in names(agent_inputs)) {
  cat(
    "\n=== RUNNING REAL",
    toupper(agent_name),
    "AGENT ===\n"
  )

  result_file <- file.path(
    test_directory,
    paste0(agent_name, "-result.rds")
  )

  if (file.exists(result_file)) {
    unlink(result_file)
  }

  command_status <- system2(
    command = "Rscript",
    args = c(
      worker_script,
      agent_name,
      agent_inputs[[agent_name]],
      result_file,
      "0.05"
    ),
    stdout = "",
    stderr = ""
  )

  if (!file.exists(result_file)) {
    validation[[agent_name]] <- data.frame(
      agent = toupper(agent_name),
      success = FALSE,
      rows = 0L,
      runtime_seconds = NA_real_,
      message = paste(
        "Worker produced no result file.",
        "Exit status:",
        command_status
      ),
      stringsAsFactors = FALSE
    )

    next
  }

  result <- readRDS(result_file)

  validation[[agent_name]] <- data.frame(
    agent = toupper(agent_name),
    success = isTRUE(result$success),
    rows = as.integer(result$rows %||% 0L),
    runtime_seconds = as.numeric(
      result$runtime_seconds %||% NA_real_
    ),
    message = as.character(
      result$message %||% result$error %||% "No message"
    ),
    stringsAsFactors = FALSE
  )
}

validation_table <- do.call(
  rbind,
  validation
)

cat("\n=== VALIDATING ALL FOUR WORKER RESULTS ===\n")

print(
  validation_table,
  row.names = FALSE
)

if (!all(validation_table$success)) {
  failed_agents <- validation_table$agent[
    !validation_table$success
  ]

  stop(
    "One or more real agents failed: ",
    paste(failed_agents, collapse = ", ")
  )
}

cat("\nREAL FOUR-AGENT WORKER TEST: PASSED\n")
