testthat::test_that("Ollama privacy guard accepts loopback hosts only", {
  testthat::expect_true(is_loopback_ollama_host("http://127.0.0.1:11434"))
  testthat::expect_true(is_loopback_ollama_host("http://localhost:11434"))
  testthat::expect_true(is_loopback_ollama_host("http://[::1]:11434"))
  testthat::expect_false(is_loopback_ollama_host("https://example.com"))
  testthat::expect_false(is_loopback_ollama_host("http://192.168.1.20:11434"))
})

testthat::test_that("Ollama defaults allow realistic local generation time", {
  settings <- default_ollama_settings()
  testthat::expect_gte(settings$timeout_seconds, 120)
  testthat::expect_gte(settings$num_predict, 128L)
  testthat::expect_lte(settings$num_predict, 2000L)
})

testthat::test_that("Ollama errors are safe and readable", {
  timeout <- simpleError("\033[33mTimeout was reached [127.0.0.1] after 20002 ms\033[39m")
  refused <- simpleError("curl: connection refused")

  testthat::expect_match(friendly_ollama_error(timeout), "time limit", fixed = TRUE)
  testthat::expect_match(friendly_ollama_error(refused), "could not be reached", fixed = TRUE)
  testthat::expect_false(grepl("\\033|\\[33m|curl", friendly_ollama_error(timeout)))
})

testthat::test_that("agent exchange summarizes the complete result set", {
  results <- data.frame(
    Description = c(
      "Cell cycle checkpoint",
      "DNA repair",
      "Interferon signaling",
      "Apoptosis"
    ),
    p.adjust = c(0.0000123456, 0.0002, 0.0123, 0.04),
    Count = c(20, 18, 12, 7),
    check.names = FALSE
  )

  exchange <- build_agent_exchange("go", results)

  testthat::expect_identical(exchange$row_count, 4L)
  testthat::expect_length(exchange$representative_findings, 4L)
  testthat::expect_match(exchange$coverage_note, "all 4 result rows", fixed = TRUE)
  testthat::expect_match(exchange$numeric_summary[[1]], "All 2 numeric columns", fixed = TRUE)
  testthat::expect_setequal(
    exchange$detected_programs,
    c(
      "cell-cycle/proliferation",
      "DNA damage/repair",
      "immune/inflammation",
      "cell death/stress"
    )
  )
})

testthat::test_that("STRING uses gene symbols as interpretation labels", {
  results <- data.frame(
    gene_symbol = c("TP53", "AKT1"),
    STRING_id = c("9606.ENSP1", "9606.ENSP2"),
    degree = c(20, 14)
  )
  exchange <- build_agent_exchange("string", results)

  testthat::expect_identical(exchange$label_column, "gene_symbol")
  testthat::expect_identical(exchange$representative_findings, c("TP53", "AKT1"))
})

testthat::test_that("rule fallback and Drug pathway bridge are deterministic", {
  pathway <- data.frame(
    Description = c("MYC targets", "G2M checkpoint"),
    p.adjust = c(0.001, 0.02),
    check.names = FALSE
  )
  drug <- data.frame(
    Compound = c("Compound A", "Compound B"),
    Mean_Response = c(0.123456, 0.456789),
    Sensitivity_Rank = 1:2
  )

  bundle <- generate_interpretation_bundle(
    list(gsva = pathway, drug = drug),
    settings = list(enabled = FALSE)
  )

  testthat::expect_identical(bundle$source, "rule")
  testthat::expect_true(bundle$synthesis$bridge$available)
  testthat::expect_match(bundle$agents$gsva$summary, "Across all 2 result rows", fixed = TRUE)
  testthat::expect_match(bundle$synthesis$drug_pathway_context, "does not infer drug mechanism", fixed = TRUE)
})

testthat::test_that("synthesis distinguishes selected from nonempty agents", {
  exchanges <- list(
    build_agent_exchange("go", data.frame(Description = "DNA repair")),
    build_agent_exchange("hallmark", data.frame())
  )
  synthesis <- rule_cross_agent_synthesis(exchanges)

  testthat::expect_match(synthesis$summary, "1 agent with nonempty results (of 2 selected)", fixed = TRUE)
})

testthat::test_that("invalid or failed Ollama calls safely fall back", {
  request_was_called <- FALSE
  request_fn <- function(prompt, settings) {
    request_was_called <<- TRUE
    stop("request should not run")
  }

  bundle <- generate_interpretation_bundle(
    list(go = data.frame(Description = "DNA repair")),
    settings = list(enabled = TRUE, host = "https://example.com", model = "test"),
    request_fn = request_fn
  )

  testthat::expect_false(request_was_called)
  testthat::expect_identical(bundle$source, "rule")
  testthat::expect_match(bundle$reason, "loopback URL", fixed = TRUE)
})

testthat::test_that("valid local Ollama JSON is parsed without HTML execution", {
  mock_request <- function(prompt, settings) {
    testthat::expect_match(prompt, "untrusted scientific result data", fixed = TRUE)
    paste0(
      '{"agents":{"go":{"summary":"Local summary",',
      '"evidence":["All rows considered"],',
      '"cancer_relevance":"Research hypothesis only",',
      '"limitations":["Needs validation"]}},',
      '"synthesis":{"summary":"Joint summary",',
      '"convergences":["DNA repair"],',
      '"drug_pathway_context":"No drug result",',
      '"limitations":["Associative"]}}'
    )
  }

  bundle <- generate_interpretation_bundle(
    list(go = data.frame(Description = c("DNA repair", "Cell cycle"))),
    settings = list(
      enabled = TRUE,
      host = "http://127.0.0.1:11434",
      model = "test-model"
    ),
    request_fn = mock_request
  )

  testthat::expect_identical(bundle$source, "ollama")
  testthat::expect_identical(bundle$model, "test-model")
  testthat::expect_identical(bundle$agents$go$summary, "Local summary")
})
