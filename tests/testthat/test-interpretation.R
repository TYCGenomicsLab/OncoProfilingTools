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
  testthat::expect_lte(settings$num_predict, 4096L)
  testthat::expect_gte(settings$num_predict, 3000L)
})

testthat::test_that("Ollama errors are safe and readable", {
  timeout <- simpleError("\033[33mTimeout was reached [127.0.0.1] after 20002 ms\033[39m")
  refused <- simpleError("curl: connection refused")

  testthat::expect_match(friendly_ollama_error(timeout), "time limit", fixed = TRUE)
  testthat::expect_match(friendly_ollama_error(refused), "could not be reached", fixed = TRUE)
  testthat::expect_false(grepl("\\033|\\[33m|curl", friendly_ollama_error(timeout)))
  testthat::expect_identical(ollama_failure_source(timeout), "ollama_timeout")
  testthat::expect_identical(ollama_failure_source(refused), "ollama_unavailable")
})

testthat::test_that("failed Ollama requests enter a terminal error state", {
  bundle <- generate_interpretation_bundle(
    list(go = data.frame(Description = "DNA repair")),
    settings = list(
      enabled = TRUE,
      host = "http://127.0.0.1:11434",
      model = "test-model"
    ),
    request_fn = function(prompt, settings) stop("connection refused")
  )

  testthat::expect_identical(bundle$source, "ollama_unavailable")
  testthat::expect_identical(bundle$model, "test-model")
  testthat::expect_match(bundle$reason, "could not be reached", fixed = TRUE)
  testthat::expect_match(bundle$agents$go$summary, "Across all 1 result rows", fixed = TRUE)
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

testthat::test_that("matrix agents expose feature-wide variation instead of early cells", {
  results <- data.frame(
    Pathway = c("Stable pathway", "Variable pathway"),
    Sample_A = c(0.1, -1),
    Sample_B = c(0.1, 0),
    Sample_C = c(0.1, 2),
    check.names = FALSE
  )
  exchange <- build_agent_exchange("gsva", results)

  testthat::expect_identical(exchange$representative_findings[[1]], "Variable pathway")
  testthat::expect_identical(exchange$feature_summaries[[1]]$feature, "Variable pathway")
  testthat::expect_identical(exchange$feature_summaries[[1]]$minimum_column, "Sample_A")
  testthat::expect_identical(exchange$feature_summaries[[1]]$maximum_column, "Sample_C")
  testthat::expect_match(build_grounded_evidence(exchange)[[1]], "across 3 measured columns", fixed = TRUE)
})

testthat::test_that("immune interpretation reports residual separately from named cells", {
  results <- data.frame(
    cell_type = c("B cell", "Macrophage M1", "T cell CD8+", "uncharacterized cell"),
    Sample_A = c(0.01, 0.04, 0.02, 0.93),
    Sample_B = c(0.03, 0.09, 0.05, 0.83),
    Sample_C = c(0.06, 0.15, 0.10, 0.69),
    check.names = FALSE
  )

  exchange <- build_agent_exchange("immune", results)
  summary <- build_matrix_observation_summary(exchange)
  evidence <- build_grounded_evidence(exchange)
  context <- rule_biological_context(exchange)

  testthat::expect_identical(exchange$row_count, 4L)
  testthat::expect_identical(exchange$interpreted_row_count, 3L)
  testthat::expect_equal(exchange$residual_compartment$median, 0.83)
  testthat::expect_false("uncharacterized cell" %in% exchange$representative_findings)
  testthat::expect_false(any(vapply(
    exchange$feature_summaries,
    function(feature) identical(feature$feature, "uncharacterized cell"),
    logical(1)
  )))
  testthat::expect_match(summary, "named immune-cell result rows", fixed = TRUE)
  testthat::expect_match(summary, "not an additional immune-cell type", fixed = TRUE)
  testthat::expect_match(evidence[[length(evidence)]], "reported separately", fixed = TRUE)
  testthat::expect_match(context, "not itself an immune-cell population", fixed = TRUE)
  testthat::expect_identical(
    sanitize_cancer_relevance(
      "Could uncharacterized cells predict treatment response?",
      agent_id = "immune"
    ),
    "Cannot be determined from the supplied result data alone."
  )
})

testthat::test_that("rich Ollama prompt requests distinct interpretation sections", {
  exchange <- build_agent_exchange("go", data.frame(Description = "DNA repair"))
  prompt <- build_ollama_prompt(list(exchange))

  testthat::expect_match(prompt, "120-200 word biological_context", fixed = TRUE)
  testthat::expect_match(prompt, "IAN-STYLE INTEGRATED REVIEW", fixed = TRUE)
  testthat::expect_match(prompt, "validation_priorities", fixed = TRUE)
  testthat::expect_match(prompt, "general biological knowledge", fixed = TRUE)
  testthat::expect_match(prompt, "Begin every biological_context exactly with", fixed = TRUE)
  testthat::expect_false(grepl("Keep each summary under 100 words", prompt, fixed = TRUE))
})

testthat::test_that("matrix summaries and biological context reject unsupported model claims", {
  results <- data.frame(
    Pathway = c("HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1"),
    Sample_A = c(-1, -0.5),
    Sample_B = c(1, 0.5),
    check.names = FALSE
  )
  exchange <- build_agent_exchange("gsva", results)
  fallback <- build_rule_interpretation_bundle(list(exchange))
  unsafe_response <- paste0(
    '{"contract_version":"', interpretation_contract_version, '","agents":{"gsva":{"summary":"The pathway is activated in 1 out of 2 samples.",',
    '"key_findings":["E2F varied"],',
    '"biological_context":"Activation suggests that cancer cells may be proliferating in these samples.",',
    '"research_hypotheses":["Test E2F"],"validation_priorities":["Compare groups"],',
    '"cancer_relevance":"Research hypothesis only.","limitations":["Associative"]}},',
    '"synthesis":{"title":"E2F profile","summary":"Single agent.","integrated_interpretation":"Single-agent profile.",',
    '"regulatory_network":"Not available.","hub_candidates":[],"convergences":[],"novelty_context":"Novelty is unverified without literature review.","next_analyses":["Compare groups"],',
    '"drug_pathway_context":"None.","limitations":["Associative"]}}'
  )
  bundle <- parse_ollama_interpretation(
    unsafe_response,
    list(exchange),
    normalise_ollama_settings(list(model = "test-model")),
    fallback
  )

  testthat::expect_identical(bundle$agents$gsva$summary, build_matrix_observation_summary(exchange))
  testthat::expect_identical(
    bundle$agents$gsva$biological_context,
    fallback$agents$gsva$biological_context
  )
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
  testthat::expect_match(bundle$agents$gsva$biological_context, "General biological context:", fixed = TRUE)
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
      '{"contract_version":"', interpretation_contract_version, '","agents":{"go":{"summary":"Across the submitted result rows, DNA repair and cell-cycle terms form the leading enrichment pattern. These named categories organize the supplied genes into related biological themes while remaining associative outputs that require an explicit comparison design and independent confirmation before direction or mechanism is assigned.",',
      '"key_findings":["DNA repair"],',
      '"biological_context":"General biological context: DNA repair terms describe systems that recognize and resolve genomic lesions, while cell-cycle terms describe replication and mitotic control. Their joint appearance can guide a testable program-level hypothesis, but enrichment does not measure pathway activity or establish a causal connection in the submitted samples.",',
      '"research_hypotheses":["Test DNA repair reproducibility"],"validation_priorities":["Independent cohort"],',
      '"cancer_relevance":"Research hypothesis only",',
      '"limitations":["Needs validation"]}},',
      '"synthesis":{"title":"DNA repair profile","summary":"Joint summary","integrated_interpretation":"Integrated DNA repair interpretation.",',
      '"regulatory_network":"No ChEA result.","hub_candidates":[],',
      '"convergences":["DNA repair"],',
      '"novelty_context":"Novelty is unverified without a literature review.","next_analyses":["Independent cohort"],',
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
  testthat::expect_identical(bundle$contract_version, interpretation_contract_version)
  testthat::expect_identical(bundle$model, "test-model")
  testthat::expect_identical(
    bundle$source_label,
    "Biological interpretation generated locally with test-model"
  )
  testthat::expect_match(bundle$agents$go$summary, "DNA repair and cell-cycle terms", fixed = TRUE)
  testthat::expect_identical(
    bundle$agents$go$biological_context,
    "General biological context: DNA repair terms describe systems that recognize and resolve genomic lesions, while cell-cycle terms describe replication and mitotic control. Their joint appearance can guide a testable program-level hypothesis, but enrichment does not measure pathway activity or establish a causal connection in the submitted samples."
  )
})
