# Background worker for one OncoProfiling analysis agent.
#
# Command:
# Rscript run_agent_worker.R AGENT INPUT_RDS RESULT_RDS PVALUE_CUTOFF

`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0) {
    fallback
  } else {
    value
  }
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript run_agent_worker.R",
      "AGENT INPUT_RDS RESULT_RDS [PVALUE_CUTOFF]"
    ),
    call. = FALSE
  )
}

agent_name <- tolower(args[[1]])
input_rds <- args[[2]]
result_rds <- args[[3]]

pvalue_cutoff <- if (length(args) >= 4) {
  suppressWarnings(as.numeric(args[[4]]))
} else {
  0.05
}

if (
  is.na(pvalue_cutoff) ||
  pvalue_cutoff <= 0 ||
  pvalue_cutoff > 1
) {
  pvalue_cutoff <- 0.05
}

all_args <- commandArgs(trailingOnly = FALSE)

file_argument <- grep(
  "^--file=",
  all_args,
  value = TRUE
)

if (length(file_argument) == 0) {
  stop(
    "Could not determine the worker script location.",
    call. = FALSE
  )
}

worker_script <- sub(
  "^--file=",
  "",
  file_argument[[1]]
)

worker_script <- normalizePath(
  worker_script,
  mustWork = TRUE
)

worker_directory <- dirname(worker_script)

project_root <- normalizePath(
  file.path(worker_directory, ".."),
  mustWork = TRUE
)

setwd(project_root)

source(
  file.path(
    project_root,
    "shiny-app",
    "run_real_agents.R"
  ),
  local = TRUE
)

if (!file.exists(input_rds)) {
  stop(
    paste(
      "Worker input file does not exist:",
      input_rds
    ),
    call. = FALSE
  )
}

input_data <- readRDS(input_rds)

started_at <- Sys.time()

cat(
  sprintf(
    "[%s] Starting %s agent.\n",
    format(started_at, "%H:%M:%S"),
    toupper(agent_name)
  )
)

worker_result <- tryCatch(
  {
    result <- run_selected_real_agent(
      agent = agent_name,
      data = input_data,
      pvalue_cutoff = pvalue_cutoff
    )

    result$started_at <- started_at
    result$finished_at <- Sys.time()

    result$runtime_seconds <- as.numeric(
      difftime(
        result$finished_at,
        result$started_at,
        units = "secs"
      )
    )

    result
  },
  error = function(error) {
    list(
      success = FALSE,
      agent = agent_name,
      rows = 0L,
      message = conditionMessage(error),
      error = conditionMessage(error),
      started_at = started_at,
      finished_at = Sys.time(),
      runtime_seconds = as.numeric(
        difftime(
          Sys.time(),
          started_at,
          units = "secs"
        )
      )
    )
  }
)


# Normalize ChEA output paths for the Shiny Results Center.
if (
  identical(tolower(agent_name), "chea") &&
  isTRUE(worker_result$success)
) {
  current_directory <- normalizePath(
    getwd(),
    mustWork = TRUE
  )

  project_root <- if (
    basename(current_directory) == "shiny-app"
  ) {
    normalizePath(
      file.path(current_directory, ".."),
      mustWork = TRUE
    )
  } else {
    current_directory
  }

  canonical_directory <- file.path(
    project_root,
    "output",
    "chea_cms4"
  )

  dir.create(
    canonical_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  canonical_csv <- file.path(
    canonical_directory,
    "chea_results.csv"
  )

  canonical_plot <- file.path(
    canonical_directory,
    "chea_tf_dotplot.png"
  )

  source_csv <- worker_result$csv %||% NA_character_
  source_plot <- worker_result$plot %||% NA_character_

  if (
    length(source_csv) == 1 &&
    !is.na(source_csv) &&
    file.exists(source_csv)
  ) {
    if (
  normalizePath(source_csv, mustWork = FALSE) !=
    normalizePath(canonical_csv, mustWork = FALSE)
) {
  file.copy(source_csv, canonical_csv, overwrite = TRUE)
}
  }

  if (
    length(source_plot) == 1 &&
    !is.na(source_plot) &&
    file.exists(source_plot)
  ) {
    if (
  normalizePath(source_plot, mustWork = FALSE) !=
    normalizePath(canonical_plot, mustWork = FALSE)
) {
  file.copy(source_plot, canonical_plot, overwrite = TRUE)
}
  }

  # Also recover paths returned inside the nested result object.
  if (
    !file.exists(canonical_csv) &&
    !is.null(worker_result$result$csv) &&
    file.exists(worker_result$result$csv)
  ) {
    file.copy(
      worker_result$result$csv,
      canonical_csv,
      overwrite = TRUE
    )
  }

  if (
    !file.exists(canonical_plot) &&
    !is.null(worker_result$result$plot) &&
    file.exists(worker_result$result$plot)
  ) {
    file.copy(
      worker_result$result$plot,
      canonical_plot,
      overwrite = TRUE
    )
  }

  worker_result$csv <- canonical_csv
  worker_result$plot <- canonical_plot

  cat(
    sprintf(
      "[%s] ChEA canonical CSV: %s | exists=%s\n",
      format(Sys.time(), "%H:%M:%S"),
      canonical_csv,
      file.exists(canonical_csv)
    )
  )

  cat(
    sprintf(
      "[%s] ChEA canonical plot: %s | exists=%s\n",
      format(Sys.time(), "%H:%M:%S"),
      canonical_plot,
      file.exists(canonical_plot)
    )
  )
}

saveRDS(
  worker_result,
  result_rds
)

if (isTRUE(worker_result$success)) {
  cat(
    sprintf(
      "[%s] %s agent completed with %s result rows.\n",
      format(Sys.time(), "%H:%M:%S"),
      toupper(agent_name),
      worker_result$rows %||% 0L
    )
  )

  quit(
    save = "no",
    status = 0L
  )
}

cat(
  sprintf(
    "[%s] %s agent failed: %s\n",
    format(Sys.time(), "%H:%M:%S"),
    toupper(agent_name),
    worker_result$message %||% "Unknown error"
  ),
  file = stderr()
)

quit(
  save = "no",
  status = 1L
)
