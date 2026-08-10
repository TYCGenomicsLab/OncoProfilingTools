testthat::test_that("every result file has a correct title and description", {
  testthat::expect_setequal(names(result_files), names(agent_titles))
  testthat::expect_setequal(names(result_files), names(agent_descriptions))
  testthat::expect_identical(agent_titles[["wikipathways"]], "WikiPathways")
  testthat::expect_identical(agent_titles[["hallmark"]], "Hallmark")
  testthat::expect_identical(agent_titles[["drug"]], "Drug Sensitivity")
  testthat::expect_identical(basename(result_files$string$plot), "string_network.png")
  testthat::expect_identical(basename(result_files$string$interactions), "string_interactions.csv")
  testthat::expect_identical(basename(result_files$immune$plot), "immune_composition_heatmap.png")
  testthat::expect_identical(basename(result_files$drug$plot), "drug_response_ranking.png")
  testthat::expect_true(all(vapply(result_files, function(files) {
    length(files$plot) == 1L && !is.na(files$plot) && nzchar(files$plot)
  }, logical(1))))
})

testthat::test_that("individual HTML reports use the correct agent title", {
  report_file <- tempfile(fileext = ".html")
  on.exit(unlink(report_file), add = TRUE)
  bundle <- generate_interpretation_bundle(list(string = data.frame()), list(enabled = FALSE))

  build_agent_html_report(report_file, "string", bundle)
  report <- paste(readLines(report_file, warn = FALSE), collapse = "\n")

  testthat::expect_match(report, "STRING Analysis Report", fixed = TRUE)
  testthat::expect_match(report, "AGENT REPORT", fixed = TRUE)
  testthat::expect_false(grepl("Cross-agent synthesis", report, fixed = TRUE))
})

testthat::test_that("display formatting does not alter source precision", {
  source_data <- data.frame(
    label = c("A", "B"),
    score = c(0.123456789, pi),
    p.adjust = c(2.916e-19, 0.0123456),
    count = c(1L, 2L)
  )
  original_score <- source_data$score
  original_p <- source_data$p.adjust
  display_data <- format_numeric_for_display(source_data)

  testthat::expect_equal(source_data$score, original_score)
  testthat::expect_equal(source_data$p.adjust, original_p)
  testthat::expect_identical(display_data$score, c("0.12", "3.14"))
  testthat::expect_identical(display_data$p.adjust, c("2.92e-19", "1.23e-02"))
  testthat::expect_identical(display_data$count, c("1.00", "2.00"))
  testthat::expect_setequal(numeric_display_columns(source_data), c("score", "p.adjust", "count"))
  testthat::expect_identical(probability_display_columns(source_data), "p.adjust")
})

testthat::test_that("Shiny UI uses namespaced small tags", {
  app_source <- paste(readLines(file.path(app_project_root, "shiny-app", "app.R"), warn = FALSE), collapse = "\n")
  results_source <- paste(readLines(file.path(app_project_root, "shiny-app", "results_helpers.R"), warn = FALSE), collapse = "\n")

  testthat::expect_false(grepl("(?<!tags\\$)small\\(", app_source, perl = TRUE))
  testthat::expect_false(grepl("(?<!tags\\$)small\\(", results_source, perl = TRUE))
})

testthat::test_that("plot token changes when file content changes", {
  plot_file <- tempfile(fileext = ".png")
  on.exit(unlink(plot_file), add = TRUE)
  writeBin(as.raw(c(1, 2, 3)), plot_file)
  first_token <- plot_cache_token(plot_file)
  writeBin(as.raw(c(3, 2, 1)), plot_file)
  second_token <- plot_cache_token(plot_file)

  testthat::expect_false(identical(first_token, second_token))
})

testthat::test_that("combined report covers all selected agents", {
  report_file <- tempfile(fileext = ".html")
  on.exit(unlink(report_file), add = TRUE)
  data_by_agent <- setNames(
    lapply(names(result_files), function(key) data.frame()),
    names(result_files)
  )
  bundle <- generate_interpretation_bundle(data_by_agent, list(enabled = FALSE))

  build_combined_html_report(
    report_file,
    interpretation_bundle = bundle,
    selected_agents = names(result_files)
  )
  report <- paste(readLines(report_file, warn = FALSE), collapse = "\n")

  testthat::expect_match(report, "WikiPathways Analysis", fixed = TRUE)
  testthat::expect_match(report, "Hallmark Analysis", fixed = TRUE)
  testthat::expect_match(report, "Immune Deconvolution Analysis", fixed = TRUE)
  testthat::expect_match(report, "Drug Sensitivity Analysis", fixed = TRUE)
  testthat::expect_match(report, "Cross-agent synthesis", fixed = TRUE)
})

testthat::test_that("Ollama loading, generating, and completion copy is explicit", {
  exchanges <- list(build_agent_exchange("go", data.frame(Description = "DNA repair")))
  loading <- build_ollama_progress_bundle(exchanges, "loading", "llama3.1:8b")
  generating <- build_ollama_progress_bundle(exchanges, "generating", "llama3.1:8b")
  completed <- generating
  completed$source <- "ollama"

  testthat::expect_identical(interpretation_source_badge(loading), "OLLAMA LOADING")
  testthat::expect_identical(
    interpretation_display_label(loading),
    "Loading llama3.1:8b locally…"
  )
  testthat::expect_identical(interpretation_source_badge(generating), "OLLAMA GENERATING")
  testthat::expect_identical(
    interpretation_display_label(generating),
    "Analyzing full result table with llama3.1:8b…"
  )
  testthat::expect_identical(interpretation_source_badge(completed), "OLLAMA COMPLETE")
  testthat::expect_identical(
    interpretation_display_label(completed),
    "Biological interpretation generated locally with llama3.1:8b"
  )
  testthat::expect_identical(interpretation_status_label(loading), "OLLAMA LOADING")
  testthat::expect_identical(interpretation_status_label(generating), "OLLAMA GENERATING")
  testthat::expect_identical(interpretation_status_label(completed), "OLLAMA COMPLETE")

  fallback <- build_rule_interpretation_bundle(exchanges)
  testthat::expect_identical(interpretation_source_badge(fallback), "SAFE FALLBACK")
  testthat::expect_identical(interpretation_status_label(fallback), "Rule summary active")
})

testthat::test_that("completed Ollama reports never retain stale progress wording", {
  report_file <- tempfile(fileext = ".html")
  on.exit(unlink(report_file), add = TRUE)
  bundle <- generate_interpretation_bundle(
    list(go = data.frame()),
    list(enabled = FALSE)
  )
  bundle$source <- "ollama"
  bundle$source_label <- "Local interpretation in progress"
  bundle$model <- "llama3.1:8b"

  build_combined_html_report(
    report_file,
    interpretation_bundle = bundle,
    selected_agents = "go"
  )
  report <- paste(readLines(report_file, warn = FALSE), collapse = "\n")

  testthat::expect_match(
    report,
    "Biological interpretation generated locally with llama3.1:8b",
    fixed = TRUE
  )
  testthat::expect_false(grepl("Local interpretation in progress", report, fixed = TRUE))
})
