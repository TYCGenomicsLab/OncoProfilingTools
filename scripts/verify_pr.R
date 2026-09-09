#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
script_path <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath("scripts/verify_pr.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Install the testthat package before running PR verification.", call. = FALSE)
}

cat("\n[1/4] Running automated tests\n")
testthat::test_dir(
  "tests/testthat",
  reporter = "summary",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE
)

cat("\n[2/4] Parsing R source files\n")
r_files <- c(
  list.files("R", "[.]R$", full.names = TRUE),
  list.files("shiny-app", "[.]R$", full.names = TRUE)
)
invisible(lapply(r_files, parse))
cat(length(r_files), "R files parsed successfully\n")

cat("\n[3/4] Checking frontend JavaScript\n")
node_status <- system2("node", c("--check", "shiny-app/www/status.js"))
if (!identical(node_status, 0L)) {
  stop("JavaScript syntax validation failed.", call. = FALSE)
}

cat("\n[4/4] Checking the Git diff\n")
git_status <- system2("git", c("diff", "--check"))
if (!identical(git_status, 0L)) {
  stop("Git diff validation failed.", call. = FALSE)
}

cat("\nAll PR verification checks passed.\n")
