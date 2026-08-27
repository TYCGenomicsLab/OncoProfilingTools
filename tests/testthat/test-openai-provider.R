openai_test_payload <- function(model = "gpt-5.6-terra") {
  payload <- list(
    contract_version = interpretation_contract_version,
    agents = list(
      go = list(
        summary = paste(
          "Across the submitted result rows, DNA repair and cell-cycle terms form a coherent enrichment pattern.",
          "These categories organize the supplied genes into related biological themes, while the result remains associative",
          "and requires an explicit comparison design plus independent validation before direction, mechanism, or disease relevance is assigned."
        ),
        key_findings = list("DNA repair"),
        biological_context = paste(
          "General biological context: DNA repair terms describe systems that recognize genomic lesions.",
          "Their appearance with cell-cycle terms supports a testable program-level hypothesis, but enrichment does not measure",
          "pathway activity or establish a causal relationship in the submitted samples."
        ),
        cancer_relevance = "Research hypothesis only",
        research_hypotheses = list("Test whether the DNA repair signal reproduces in an independent cohort."),
        validation_priorities = list("Independent cohort"),
        limitations = list("Enrichment is associative and database dependent.")
      )
    ),
    synthesis = list(
      title = "DNA repair profile",
      summary = "A result-grounded DNA repair summary.",
      integrated_interpretation = "The deterministic findings support a bounded DNA repair hypothesis.",
      regulatory_network = "No regulatory result was supplied.",
      hub_candidates = list(),
      convergences = list("DNA repair"),
      novelty_context = "Novelty is unverified without a literature review.",
      next_analyses = list("Independent cohort"),
      drug_pathway_context = "No drug result was supplied.",
      limitations = list("Associative result")
    )
  )
  jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null")
}

with_test_openai_key <- function(code) {
  old <- Sys.getenv("OPENAI_API_KEY", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("OPENAI_API_KEY") else Sys.setenv(OPENAI_API_KEY = old)
  }, add = TRUE)
  Sys.setenv(OPENAI_API_KEY = "sk-test-never-serialize-this-value")
  force(code)
}

testthat::test_that("OpenAI settings require an environment key and explicit data consent", {
  old <- Sys.getenv("OPENAI_API_KEY", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("OPENAI_API_KEY") else Sys.setenv(OPENAI_API_KEY = old)
  }, add = TRUE)

  Sys.unsetenv("OPENAI_API_KEY")
  testthat::expect_match(
    validate_openai_settings(list(provider = "openai", openai_data_consent = TRUE)),
    "not configured",
    fixed = TRUE
  )

  Sys.setenv(OPENAI_API_KEY = "sk-test-never-serialize-this-value")
  testthat::expect_match(
    validate_openai_settings(list(provider = "openai", openai_data_consent = FALSE)),
    "Confirm",
    fixed = TRUE
  )
  testthat::expect_null(validate_openai_settings(list(provider = "openai", openai_data_consent = TRUE)))

  settings <- normalise_interpretation_settings(list(provider = "openai", openai_data_consent = TRUE))
  captured <- paste(capture.output(str(settings)), collapse = " ")
  testthat::expect_false(grepl("sk-test-never-serialize-this-value", captured, fixed = TRUE))
  testthat::expect_true(settings$openai_key_available)
})

testthat::test_that("Responses API text extraction handles canonical output content", {
  payload <- list(output = list(list(content = list(
    list(type = "reasoning", text = "hidden"),
    list(type = "output_text", text = "grounded response")
  ))))
  testthat::expect_identical(extract_openai_response_text(payload), "grounded response")
})

testthat::test_that("OpenAI strict schema declares only the selected scientific agents", {
  schema <- openai_interpretation_schema(c("go", "string"))
  agent_schema <- schema$properties$agents

  testthat::expect_false(agent_schema$additionalProperties)
  testthat::expect_setequal(names(agent_schema$properties), c("go", "string"))
  testthat::expect_setequal(agent_schema$required, c("go", "string"))
  testthat::expect_true(all(vapply(
    agent_schema$properties,
    function(value) identical(value$additionalProperties, FALSE),
    logical(1)
  )))
  testthat::expect_error(openai_interpretation_schema(character()), "at least one")
})

testthat::test_that("OpenAI interpretation uses the same grounding contract without a network call", {
  with_test_openai_key({
    observed_prompt <- NULL
    mock_openai <- function(prompt, settings) {
      observed_prompt <<- prompt
      value <- openai_test_payload()
      attr(value, "openai_metadata") <- list(
        request_id = "resp_test_123",
        model = settings$openai_model,
        usage = list(input_tokens = 800, output_tokens = 400, total_tokens = 1200)
      )
      value
    }
    bundle <- generate_openai_interpretation_bundle(
      list(go = data.frame(Description = c("DNA repair", "Cell cycle"))),
      settings = list(provider = "openai", openai_data_consent = TRUE),
      request_fn = mock_openai
    )

    testthat::expect_identical(bundle$source, "openai")
    testthat::expect_identical(bundle$provider, "openai")
    testthat::expect_identical(bundle$request_id, "resp_test_123")
    testthat::expect_identical(bundle$usage$total_tokens, 1200)
    testthat::expect_true(is.finite(bundle$estimated_cost_usd))
    testthat::expect_match(observed_prompt, "untrusted scientific result data", fixed = TRUE)
    testthat::expect_match(observed_prompt, "DATA_START", fixed = TRUE)
    testthat::expect_identical(attr(observed_prompt, "openai_agent_ids"), "go")
    testthat::expect_identical(bundle$agents$go$key_findings, c("DNA repair", "Cell cycle"))
  })
})

testthat::test_that("Compare mode preserves both providers and selects OpenAI when both complete", {
  with_test_openai_key({
    local_mock <- function(prompt, settings) openai_test_payload()
    remote_mock <- function(prompt, settings) {
      value <- openai_test_payload()
      attr(value, "openai_metadata") <- list(
        request_id = "resp_compare_456",
        model = settings$openai_model,
        usage = list(input_tokens = 700, output_tokens = 300, total_tokens = 1000)
      )
      value
    }
    bundle <- generate_provider_interpretation_bundle(
      list(go = data.frame(Description = c("DNA repair", "Cell cycle"))),
      settings = list(
        provider = "compare",
        enabled = TRUE,
        host = "http://127.0.0.1:11434",
        model = "llama3.1:8b",
        openai_data_consent = TRUE
      ),
      ollama_request_fn = local_mock,
      openai_request_fn = remote_mock
    )

    testthat::expect_identical(bundle$source, "comparison")
    testthat::expect_identical(bundle$comparison$primary_provider, "openai")
    testthat::expect_identical(bundle$comparison$ollama$source, "ollama")
    testthat::expect_identical(bundle$comparison$openai$source, "openai")
    comparison_html <- report_provider_comparison_html(bundle, "go")
    testthat::expect_match(comparison_html, "Ollama vs OpenAI Premium", fixed = TRUE)
    testthat::expect_match(comparison_html, "store=false", fixed = TRUE)
    testthat::expect_false(grepl("sk-test-never-serialize-this-value", comparison_html, fixed = TRUE))
  })
})

testthat::test_that("OpenAI failures terminate honestly with deterministic observations", {
  with_test_openai_key({
    bundle <- generate_openai_interpretation_bundle(
      list(go = data.frame(Description = "DNA repair")),
      settings = list(provider = "openai", openai_data_consent = TRUE),
      request_fn = function(prompt, settings) stop("401 invalid token")
    )
    testthat::expect_true(bundle$source %in% c("openai_error", "openai_unavailable"))
    testthat::expect_true(length(bundle$agents$go$observed_results) > 0L)
    testthat::expect_false(grepl("sk-test-never-serialize-this-value", bundle$reason, fixed = TRUE))
  })
})

testthat::test_that("OpenAI HTTP failures retain safe actionable diagnostics", {
  with_test_openai_key({
    bundle <- generate_openai_interpretation_bundle(
      list(go = data.frame(Description = "DNA repair")),
      settings = list(provider = "openai", openai_data_consent = TRUE),
      request_fn = function(prompt, settings) {
        stop(new_openai_api_error(
          "Invalid schema containing sk-secret-should-not-escape",
          http_status = 400L,
          request_id = "req_schema_123"
        ))
      }
    )

    testthat::expect_identical(bundle$source, "openai_error")
    testthat::expect_identical(bundle$http_status, 400L)
    testthat::expect_identical(bundle$request_id, "req_schema_123")
    testthat::expect_match(bundle$reason, "structured-output request", fixed = TRUE)
    testthat::expect_match(bundle$provider_error, "[REDACTED API KEY]", fixed = TRUE)
    testthat::expect_false(grepl("sk-secret-should-not-escape", bundle$provider_error, fixed = TRUE))
  })
})

testthat::test_that("failed provider comparisons do not present computed prose as OpenAI output", {
  with_test_openai_key({
    bundle <- generate_provider_interpretation_bundle(
      list(go = data.frame(Description = c("DNA repair", "Cell cycle"))),
      settings = list(
        provider = "compare",
        enabled = TRUE,
        host = "http://127.0.0.1:11434",
        model = "llama3.1:8b",
        openai_data_consent = TRUE
      ),
      ollama_request_fn = function(prompt, settings) openai_test_payload(),
      openai_request_fn = function(prompt, settings) {
        stop(new_openai_api_error(
          "Invalid structured-output schema",
          http_status = 400L,
          request_id = "req_failed_comparison"
        ))
      }
    )

    html <- report_provider_comparison_html(bundle, "go")
    testthat::expect_match(html, "OpenAI did not generate model-authored prose", fixed = TRUE)
    testthat::expect_match(html, "req_failed_comparison", fixed = TRUE)
    testthat::expect_match(html, "HTTP status", fixed = TRUE)
  })
})
