# Real background-pipeline helpers for the Shiny application.

real_pipeline_project_root <- local({
  current_directory <- normalizePath(
    getwd(),
    mustWork = TRUE
  )

  if (
    basename(current_directory) == "shiny-app"
  ) {
    normalizePath(
      file.path(current_directory, ".."),
      mustWork = TRUE
    )
  } else if (
    dir.exists(
      file.path(current_directory, "shiny-app")
    )
  ) {
    current_directory
  } else {
    stop(
      "Could not locate the OncoProfilingTools project root.",
      call. = FALSE
    )
  }
})

real_pipeline_runtime_dir <- file.path(
  real_pipeline_project_root,
  "shiny-app",
  "runtime"
)

dir.create(
  real_pipeline_runtime_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

real_pipeline_worker_script <- file.path(
  real_pipeline_project_root,
  "shiny-app",
  "run_agent_worker.R"
)

create_real_analysis_run <- function(data) {
  run_id <- paste0(
    format(Sys.time(), "%Y%m%d-%H%M%S"),
    "-",
    sample.int(999999, 1)
  )

  run_directory <- file.path(
    real_pipeline_runtime_dir,
    run_id
  )

  dir.create(
    run_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  input_file <- file.path(
    run_directory,
    "input.rds"
  )

  saveRDS(
    data,
    input_file
  )

  list(
    id = run_id,
    directory = run_directory,
    input_file = input_file
  )
}

start_real_agent_worker <- function(
  agent,
  analysis_run,
  pvalue_cutoff = 0.05
) {
  agent <- tolower(agent)

  result_file <- file.path(
    analysis_run$directory,
    paste0(agent, "-result.rds")
  )

  log_file <- file.path(
    analysis_run$directory,
    paste0(agent, ".log")
  )

  process <- processx::process$new(
    command = file.path(
      R.home("bin"),
      "Rscript"
    ),
    args = c(
      real_pipeline_worker_script,
      agent,
      analysis_run$input_file,
      result_file,
      as.character(pvalue_cutoff)
    ),
    stdout = log_file,
    stderr = "2>&1",
    cleanup = TRUE
  )

  list(
    agent = agent,
    process = process,
    result_file = result_file,
    log_file = log_file
  )
}

read_real_worker_result <- function(worker) {
  if (!file.exists(worker$result_file)) {
    return(NULL)
  }

  tryCatch(
    readRDS(worker$result_file),
    error = function(error) {
      NULL
    }
  )
}

read_real_worker_log <- function(worker) {
  if (!file.exists(worker$log_file)) {
    return(character())
  }

  readLines(
    worker$log_file,
    warn = FALSE
  )
}
