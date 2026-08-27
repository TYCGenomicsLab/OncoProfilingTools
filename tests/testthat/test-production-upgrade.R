testthat::test_that("1 headerless Ensembl text is imported as data", {
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("ENSG00000141510", "ENSG00000146648"), path)
  data <- read_analysis_dataset(path, basename(path))
  testthat::expect_identical(names(data), "gene_id")
  testthat::expect_false(attr(data, "import_metadata")$header_detected)
  testthat::expect_equal(nrow(data), 2L)
})

testthat::test_that("2 named DEG headers are retained", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("Geneid,FDR,logFC", "TP53,0.01,2"), path)
  data <- read_analysis_dataset(path, basename(path))
  testthat::expect_true(attr(data, "import_metadata")$header_detected)
  testthat::expect_identical(names(data), c("Geneid", "FDR", "logFC"))
})

testthat::test_that("3 Ensembl versions are normalized before mapping", {
  mapping <- map_gene_identifiers(c("ENSG00000141510.18", "ENSG00000146648.22"))
  testthat::expect_identical(mapping$metadata$identifier_type, "ensembl")
  testthat::expect_true(all(c("TP53", "EGFR") %in% mapping$symbols))
})

testthat::test_that("4 unmapped identifiers are counted honestly", {
  mapping <- map_gene_identifiers(c("ENSG00000141510", "ENSG99999999999"), "ensembl")
  testthat::expect_equal(mapping$metadata$input_count, 2L)
  testthat::expect_equal(mapping$metadata$mapped_count, 1L)
  testthat::expect_equal(mapping$metadata$unmapped_count, 1L)
  testthat::expect_identical(mapping$metadata$unmapped_examples, "ENSG99999999999")
})

testthat::test_that("5 IAN-style FDR and logFC filters are applied", {
  data <- data.frame(
    Geneid = c("TP53", "EGFR", "BRCA1"),
    FDR = c(0.01, 0.2, 0.03),
    logFC = c(2, 3, 0.2)
  )
  prepared <- prepare_gene_input(data, pvalue_cutoff = 0.05, effect_cutoff = 1)
  testthat::expect_identical(prepared$genes, c("TP53", "BRCA1"))
  testthat::expect_match(prepared$selection_note, "FDR ≤ 0.05", fixed = TRUE)
})

testthat::test_that("6 selected DEG rows are capped reproducibly", {
  data <- data.frame(
    Geneid = paste0("GENE", 1:10),
    FDR = seq(0.001, 0.01, length.out = 10),
    logFC = seq(2, 3, length.out = 10)
  )
  prepared <- prepare_gene_input(data, max_genes = 4L)
  testthat::expect_equal(length(prepared$genes), 4L)
  testthat::expect_match(prepared$selection_note, "capped", fixed = TRUE)
})

testthat::test_that("7 near-whole-genome unranked lists are rejected", {
  data <- data.frame(gene = paste0("GENE", seq_len(5001)))
  testthat::expect_error(prepare_gene_input(data), "not a meaningful over-representation input", fixed = TRUE)
})

testthat::test_that("8 symbol inputs keep one-to-one mapping provenance", {
  prepared <- prepare_gene_input(data.frame(SYMBOL = c("tp53", "EGFR", "PTEN")))
  testthat::expect_identical(prepared$mapping_metadata$identifier_type, "symbol")
  testthat::expect_equal(prepared$mapping_metadata$mapping_rate, 1)
  testthat::expect_identical(prepared$genes, c("TP53", "EGFR", "PTEN"))
})

testthat::test_that("9 Ensembl expression matrices map to symbol row names", {
  ids <- c(
    "ENSG00000141510", "ENSG00000146648", "ENSG00000139618", "ENSG00000012048",
    "ENSG00000171862", "ENSG00000157764", "ENSG00000133703", "ENSG00000141736",
    "ENSG00000142208", "ENSG00000121879", "ENSG00000136997", "ENSG00000105221"
  )
  data <- data.frame(gene_id = ids, matrix(seq_len(120), nrow = 12), check.names = FALSE)
  names(data)[-1] <- paste0("sample_", 1:10)
  matrix <- prepare_expression_matrix(data)
  testthat::expect_true("TP53" %in% rownames(matrix))
  testthat::expect_equal(ncol(matrix), 10L)

  wide <- data.frame(sample_id = paste0("S", 1:10), matrix(seq_len(120), nrow = 10), check.names = FALSE)
  names(wide)[-1] <- paste0("GENE", 1:12, " (", 101:112, ")")
  wide_matrix <- prepare_expression_matrix(wide)
  testthat::expect_equal(dim(wide_matrix), c(12L, 10L))
  testthat::expect_identical(colnames(wide_matrix), paste0("S", 1:10))
  testthat::expect_true("GENE1" %in% rownames(wide_matrix))
})

testthat::test_that("10 expression matrices require multiple samples", {
  data <- data.frame(gene = paste0("GENE", 1:12), sample = seq_len(12))
  testthat::expect_error(prepare_expression_matrix(data), "at least two numeric sample columns", fixed = TRUE)
})


testthat::test_that("11 landing page exposes two separate workflows", {
  source_text <- paste(readLines(file.path(app_project_root, "shiny-app", "workflow_ui.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(source_text, "Biomarker Discovery", fixed = TRUE)
  testthat::expect_match(source_text, "Drug Sensitivity", fixed = TRUE)
  testthat::expect_match(source_text, "nav_biomarker", fixed = TRUE)
  testthat::expect_match(source_text, "nav_drug", fixed = TRUE)
  testthat::expect_false(grepl("(?<!tags\\$)small\\(", source_text, perl = TRUE))
})

testthat::test_that("12 biomarker and drug agent keys are disjoint", {
  source_text <- paste(readLines(file.path(app_project_root, "shiny-app", "workflow_ui.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(source_text, 'setdiff\\(names\\(module_meta\\), "drug"\\)')
  testthat::expect_match(source_text, 'drug_agent_keys <- "drug"', fixed = TRUE)
})

testthat::test_that("13 oversized duplicate summary UI is absent from the route shell", {
  source_text <- paste(readLines(file.path(app_project_root, "shiny-app", "workflow_ui.R"), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("final_summary|progress-section|summary-modules", source_text))
  testthat::expect_match(source_text, "run_status_bar", fixed = TRUE)
})

testthat::test_that("14 pastel layer defines a 60-30-10 style palette", {
  css <- paste(readLines(file.path(app_project_root, "shiny-app", "www", "pastel.css"), warn = FALSE), collapse = "\n")
  testthat::expect_match(css, "#f7f5ef", fixed = TRUE)
  testthat::expect_match(css, "#dcebdd", fixed = TRUE)
  testthat::expect_match(css, "#bd6049", fixed = TRUE)
  testthat::expect_match(css, "prefers-reduced-motion", fixed = TRUE)
})

testthat::test_that("15 results center labels workflows independently", {
  source_text <- paste(readLines(file.path(app_project_root, "shiny-app", "results_helpers.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(source_text, "Biomarker Results Center", fixed = TRUE)
  testthat::expect_match(source_text, "Drug Sensitivity Results Center", fixed = TRUE)
  testthat::expect_match(source_text, "length(selected_agents) >= 1L", fixed = TRUE)
})

testthat::test_that("16 interpretation schema requires the versioned contract", {
  schema <- ollama_interpretation_schema()
  testthat::expect_true("contract_version" %in% schema$required)
  testthat::expect_identical(schema$properties$contract_version$enum[[1L]], interpretation_contract_version)
})

testthat::test_that("17 missing interpretation contract is rejected", {
  exchange <- build_agent_exchange("go", data.frame(Description = "DNA repair"))
  fallback <- build_rule_interpretation_bundle(list(exchange))
  response <- '{"agents":{},"synthesis":{}}'
  testthat::expect_error(
    parse_ollama_interpretation(response, list(exchange), normalise_ollama_settings(), fallback),
    "contract mismatch",
    fixed = TRUE
  )
})

testthat::test_that("18 deterministic observations are present in every fallback entry", {
  bundle <- build_rule_interpretation_bundle(list(build_agent_exchange("go", data.frame(Description = "DNA repair"))))
  testthat::expect_identical(bundle$agents$go$observed_results, bundle$agents$go$evidence)
  testthat::expect_true(length(bundle$agents$go$observed_results) > 0L)
  testthat::expect_identical(
    sanitize_cancer_relevance("These findings may have implications for cancer treatment."),
    "Cannot be determined from the supplied result data alone."
  )
  drug_exchange <- build_agent_exchange("drug", data.frame(Compound = "Example", AUC = 0.4))
  testthat::expect_identical(
    sanitize_biological_context(
      "drug",
      "General biological context: Example is a proteasome inhibitor with a known target.",
      rule_biological_context(drug_exchange)
    ),
    rule_biological_context(drug_exchange)
  )
})

testthat::test_that("19 production report records provenance mapping and manifest", {
  report <- tempfile(fileext = ".html")
  on.exit(unlink(report), add = TRUE)
  bundle <- generate_interpretation_bundle(list(go = data.frame()), list(enabled = FALSE))
  context <- list(
    workflow = "biomarker",
    input = list(name = "deg.txt", checksum_md5 = "abc123", rows = 180L, columns = 1L),
    mapping = list(identifier_type = "ensembl", input_count = 180L, mapped_count = 137L, output_symbol_count = 137L, unmapped_count = 43L, mapping_rate = 137 / 180, duplicate_mappings_removed = 0L)
  )
  build_combined_html_report(report, bundle, "go", run_context = context)
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  testthat::expect_match(html, "deg.txt", fixed = TRUE)
  testthat::expect_match(html, "76.1%", fixed = TRUE)
  testthat::expect_match(html, "ARTIFACT MANIFEST", fixed = TRUE)
  testthat::expect_lt(regexpr("Responsible interpretation", html, fixed = TRUE)[[1L]], regexpr("IAN-STYLE INTEGRATED REVIEW", html, fixed = TRUE)[[1L]])
  testthat::expect_gt(regexpr("Gene mapping summary", html, fixed = TRUE)[[1L]], regexpr("id='methods'", html, fixed = TRUE)[[1L]])
})

testthat::test_that("20 progress reports are forced to a terminal export state", {
  report <- tempfile(fileext = ".html")
  on.exit(unlink(report), add = TRUE)
  exchanges <- list(build_agent_exchange("go", data.frame(Description = "DNA repair")))
  progress <- build_ollama_progress_bundle(exchanges, "generating", "llama3.1:8b")
  build_combined_html_report(report, progress, "go")
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  testthat::expect_match(html, "Computed scientific summary", fixed = TRUE)
  testthat::expect_false(grepl("OLLAMA GENERATING|Analyzing full result table", html))
})

testthat::test_that("21 interactive visual summaries produce professional bar and 3D metrics", {
  data <- data.frame(
    Description = c("DNA repair", "Cell cycle", "Apoptosis"),
    p.adjust = c(1e-6, 1e-4, 0.01),
    Count = c(18, 12, 7),
    FoldEnrichment = c(3.2, 2.4, 1.8)
  )
  summary <- result_visual_summary(data, "go")
  testthat::expect_equal(nrow(summary$data), 3L)
  testthat::expect_match(summary$metric, "−log10", fixed = TRUE)
  bar_plot <- professional_bar_plot(summary, "GO")
  testthat::expect_s3_class(bar_plot, "plotly")
  built <- plotly::plotly_build(bar_plot)
  testthat::expect_true(length(built$x$data) > 0L)
  testthat::expect_false(grepl("Tealgrn", paste(unlist(built$x$data), collapse = " "), fixed = TRUE))
  testthat::expect_null(evidence_3d_plot(summary, "GO"))
})

testthat::test_that("21c STRING 3D network uses retrieved edges", {
  hubs <- data.frame(gene_symbol = c("EGFR", "FN1", "CD44", "CDH1"), degree = c(206, 149, 145, 143))
  interactions <- data.frame(
    from_name = c("EGFR", "EGFR", "FN1", "CD44"),
    to_name = c("FN1", "CD44", "CDH1", "CDH1"),
    combined_score = c(944, 998, 873, 927)
  )
  plot <- evidence_3d_plot(result_visual_summary(hubs, "string"), "STRING", interactions)
  built <- plotly::plotly_build(plot)
  testthat::expect_true(length(built$x$data) >= 2L)
  testthat::expect_identical(built$x$data[[1L]]$mode, "lines")
  testthat::expect_match(built$x$data[[2L]]$mode, "markers+text", fixed = TRUE)
  testthat::expect_setequal(built$x$data[[2L]]$text, hubs$gene_symbol)
})

testthat::test_that("21b agent evidence cannot cross-contaminate ChEA interpretation", {
  chea <- build_agent_exchange(
    "chea",
    data.frame(
      Term = c("SOX2 20726797 ChIP-Seq SW620 Human", "SUZ12 20007587 ChIP-Seq Mouse", "NR3C1 27076634 ChIP-Seq BEAS2B Human"),
      P.value = c(1e-40, 4.6e-25, 3.5e-22),
      Combined.Score = c(1507, 296.7, 139),
      check.names = FALSE
    )
  )
  fallback <- build_rule_interpretation_bundle(list(chea))
  response <- jsonlite::toJSON(list(
    contract_version = interpretation_contract_version,
    agents = list(chea = list(
      summary = "Transcriptional regulation analysis results",
      key_findings = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_KRAS_SIGNALING_UP"),
      biological_context = "General biological context: Hallmark estrogen and KRAS programs describe broad cancer biology.",
      cancer_relevance = "These findings are relevant to breast cancer, lung cancer, and colon cancer.",
      research_hypotheses = c("Test HALLMARK_KRAS_SIGNALING_UP", "Validate estrogen-response signaling"),
      validation_priorities = c("Repeat the analysis", "Use an orthogonal assay", "Confirm mapping"),
      limitations = c("Associative result")
    )),
    synthesis = fallback$synthesis
  ), auto_unbox = TRUE)
  parsed <- parse_ollama_interpretation(
    response, list(chea), normalise_ollama_settings(list(model = "test-model")), fallback
  )
  entry <- parsed$agents$chea
  testthat::expect_identical(entry$key_findings, utils::head(chea$representative_findings, 6L))
  testthat::expect_false(any(grepl("HALLMARK", entry$key_findings, fixed = TRUE)))
  testthat::expect_identical(entry$summary, fallback$agents$chea$summary)
  testthat::expect_identical(entry$biological_context, fallback$agents$chea$biological_context)
  testthat::expect_identical(entry$research_hypotheses, fallback$agents$chea$research_hypotheses)
  testthat::expect_identical(entry$cancer_relevance, "Cannot be determined from the supplied result data alone.")
  testthat::expect_match(build_grounded_evidence(chea)[[1L]], "1e-40", fixed = TRUE)
})

testthat::test_that("22 user-facing source copy never says safe fallback", {
  results_text <- paste(readLines(file.path(app_project_root, "shiny-app", "results_helpers.R"), warn = FALSE), collapse = "\n")
  workflow_text <- paste(readLines(file.path(app_project_root, "shiny-app", "workflow_ui.R"), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("SAFE FALLBACK|Rule-based fallback", results_text, ignore.case = TRUE))
  testthat::expect_false(grepl("SAFE FALLBACK|Rule-based fallback", workflow_text, ignore.case = TRUE))
  testthat::expect_match(results_text, "No provisional computed narrative is shown", fixed = TRUE)
})

testthat::test_that("23 detailed report embeds IAN review and interactive evidence graphs", {
  report <- tempfile(fileext = ".html")
  on.exit(unlink(report), add = TRUE)
  original <- result_files$go$csv
  fixture <- tempfile(fileext = ".csv")
  readr::write_csv(data.frame(Description = c("DNA repair", "Cell cycle", "Apoptosis"), p.adjust = c(1e-6, 1e-4, 0.01), Count = c(18, 12, 7), FoldEnrichment = c(3.2, 2.4, 1.8)), fixture)
  result_files$go$csv <<- fixture
  on.exit(result_files$go$csv <<- original, add = TRUE)
  bundle <- generate_interpretation_bundle(list(go = safe_result_csv(fixture)), list(enabled = FALSE))
  build_combined_html_report(report, bundle, "go")
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  testthat::expect_match(html, "IAN-STYLE INTEGRATED REVIEW", fixed = TRUE)
  testthat::expect_match(html, "Plotly.newPlot", fixed = TRUE)
  testthat::expect_false(grepl("Recommended validation priorities", html, fixed = TRUE))
  testthat::expect_match(html, "Biological and cellular interpretation", fixed = TRUE)
  testthat::expect_match(html, "One consolidated result-grounded hypothesis", fixed = TRUE)
  testthat::expect_false(grepl("plotly-3d-", html, fixed = TRUE))
  testthat::expect_false(grepl("Full-precision result preview", html, fixed = TRUE))
  testthat::expect_false(grepl("Result-grounded research hypotheses", html, fixed = TRUE))
  testthat::expect_false(grepl("<h3>Recommended next analyses", html, fixed = TRUE))
  testthat::expect_false(grepl("<h3>Interpretation limits", html, fixed = TRUE))
})

testthat::test_that("23b report keeps a documented 3D network for STRING only", {
  report <- tempfile(fileext = ".html")
  nodes_file <- tempfile(fileext = ".csv")
  interactions_file <- tempfile(fileext = ".csv")
  on.exit(unlink(c(report, nodes_file, interactions_file)), add = TRUE)
  nodes <- data.frame(
    STRING_id = paste0("9606.ENSP", 1:4),
    gene_symbol = c("EGFR", "FN1", "CD44", "CDH1"),
    degree = c(206, 149, 145, 143)
  )
  interactions <- data.frame(
    from_name = c("EGFR", "EGFR", "FN1", "CD44"),
    to_name = c("FN1", "CD44", "CDH1", "CDH1"),
    combined_score = c(944, 998, 873, 927)
  )
  readr::write_csv(nodes, nodes_file)
  readr::write_csv(interactions, interactions_file)
  original_nodes <- result_files$string$csv
  original_interactions <- result_files$string$interactions
  result_files$string$csv <<- nodes_file
  result_files$string$interactions <<- interactions_file
  on.exit({
    result_files$string$csv <<- original_nodes
    result_files$string$interactions <<- original_interactions
  }, add = TRUE)
  bundle <- generate_interpretation_bundle(list(string = nodes), list(enabled = FALSE))
  build_combined_html_report(report, bundle, "string")
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  testthat::expect_match(html, "plotly-3d-string", fixed = TRUE)
  testthat::expect_match(html, "Connected 3D STRING interaction network", fixed = TRUE)
  testthat::expect_match(html, "NETWORK GUIDE", fixed = TRUE)
  testthat::expect_match(html, "eligible retrieved edges", fixed = TRUE)
  testthat::expect_match(html, "Highest-degree hubs", fixed = TRUE)
  testthat::expect_match(html, "Drag</strong> to rotate", fixed = TRUE)
  testthat::expect_match(html, "connectivity prioritizes candidates", fixed = TRUE)
})

testthat::test_that("24 results navigation exposes a prominent return to the main menu", {
  workflow_text <- paste(readLines(file.path(app_project_root, "shiny-app", "workflow_ui.R"), warn = FALSE), collapse = "\n")
  results_text <- paste(readLines(file.path(app_project_root, "shiny-app", "results_helpers.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(workflow_text, "Back to main menu", fixed = TRUE)
  testthat::expect_match(results_text, "results_back_home", fixed = TRUE)
})
