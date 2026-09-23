app_project_root <- normalizePath(
  testthat::test_path("..", ".."),
  mustWork = TRUE
)

old_test_working_directory <- getwd()
setwd(app_project_root)
source(file.path("shiny-app", "interpretation.R"), local = TRUE)
source(file.path("shiny-app", "interactive_visuals.R"), local = TRUE)
source(file.path("shiny-app", "results_helpers.R"), local = TRUE)
source(file.path("shiny-app", "production_reporting.R"), local = TRUE)
source(file.path("shiny-app", "run_real_agents.R"), local = TRUE)
setwd(old_test_working_directory)
