library(testthat)
library(OncoProfilingTools)

# The Shiny application is intentionally not installed as part of the R
# package. Run these source-level tests from the repository with
# testthat::test_dir("tests/testthat"). During an installed-package check the
# source tree is unavailable, so the app tests are reported by the dedicated
# repository command instead.
if (dir.exists(file.path("..", "shiny-app"))) {
  test_check("OncoProfilingTools")
} else {
  message("Shiny source tests: use testthat::test_dir('tests/testthat') from the repository root.")
}
