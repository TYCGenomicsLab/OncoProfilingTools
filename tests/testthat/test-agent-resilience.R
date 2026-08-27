testthat::test_that("DEG tables are converted to bounded analysis gene sets", {
  data <- data.frame(
    hugo_symbol = paste0("GENE", seq_len(6000)),
    entrez_id = seq_len(6000),
    diff_case_minus_control = seq(-3, 3, length.out = 6000),
    stringsAsFactors = FALSE
  )

  prepared <- prepare_gene_input(data, max_genes = 1000L)

  testthat::expect_length(prepared$genes, 1000L)
  testthat::expect_equal(nrow(prepared$data), 1000L)
  testthat::expect_identical(prepared$effect_column, "diff_case_minus_control")
  testthat::expect_match(prepared$selection_note, "1,000 of 6,000", fixed = TRUE)
  testthat::expect_match(prepared$selection_note, "positive and negative effects are combined", fixed = TRUE)
})

testthat::test_that("very large unranked lists fail with actionable guidance", {
  data <- data.frame(gene_symbol = paste0("GENE", seq_len(5001)))
  testthat::expect_error(
    prepare_gene_input(data),
    "near-whole-genome list",
    fixed = TRUE
  )
})

testthat::test_that("ChEA retries expired responses and accepts only the real schema", {
  testthat::skip_if_not_installed("enrichR")
  attempts <- 0L
  mock_request <- function(genes, database) {
    attempts <<- attempts + 1L
    if (attempts == 1L) {
      return(stats::setNames(list(data.frame(X. = c('"expired": true', "}"))), database))
    }
    stats::setNames(
      list(data.frame(
        Term = "TP53 ChIP-seq",
        Adjusted.P.value = 0.001,
        Combined.Score = 42,
        check.names = FALSE
      )),
      database
    )
  }

  result <- run_chea_agent(
    c("TP53", "CDKN1A"),
    request_fn = mock_request,
    sleep_fn = function(seconds) NULL,
    retry_delay_seconds = 0
  )

  testthat::expect_identical(attempts, 2L)
  testthat::expect_identical(result$results$Term, "TP53 ChIP-seq")
})

testthat::test_that("ChEA never saves an expired payload as science", {
  testthat::skip_if_not_installed("enrichR")
  expired_request <- function(genes, database) {
    stats::setNames(list(data.frame(X. = c('"expired": true', "}"))), database)
  }

  testthat::expect_error(
    run_chea_agent(
      c("TP53", "CDKN1A"),
      max_attempts = 2L,
      request_fn = expired_request,
      sleep_fn = function(seconds) NULL,
      retry_delay_seconds = 0
    ),
    "expired API data was not saved",
    fixed = TRUE
  )
})

testthat::test_that("KEGG retries transient REST failures with bounded backoff", {
  attempts <- 0L
  delays <- numeric()
  mock_enrichment <- function(...) {
    attempts <<- attempts + 1L
    if (attempts < 3L) stop("cannot read from connection", call. = FALSE)
    structure(list(ok = TRUE), class = "mock_kegg_result")
  }

  result <- run_kegg_enrichment_with_retry(
    gene = c("1956", "7157"),
    max_attempts = 3L,
    retry_delay_seconds = 1,
    request_timeout_seconds = 90,
    enrich_fn = mock_enrichment,
    sleep_fn = function(seconds) delays <<- c(delays, seconds)
  )

  testthat::expect_s3_class(result, "mock_kegg_result")
  testthat::expect_identical(attempts, 3L)
  testthat::expect_identical(delays, c(1, 2))
})

testthat::test_that("KEGG does not retry permanent input errors", {
  attempts <- 0L
  testthat::expect_error(
    run_kegg_enrichment_with_retry(
      gene = "1956",
      max_attempts = 3L,
      retry_delay_seconds = 0,
      enrich_fn = function(...) {
        attempts <<- attempts + 1L
        stop("invalid organism code", call. = FALSE)
      },
      sleep_fn = function(seconds) NULL
    ),
    "after 1 attempt",
    fixed = TRUE
  )
  testthat::expect_identical(attempts, 1L)
})

testthat::test_that("STRING creates a visualization from network results", {
  plot_file <- tempfile(fileext = ".png")
  on.exit(unlink(plot_file), add = TRUE)

  result <- list(
    nodes = data.frame(
      STRING_id = paste0("9606.ENSP", 1:4),
      gene_symbol = c("TP53", "AKT1", "MYC", "EGFR"),
      degree = c(3, 2, 2, 1),
      stringsAsFactors = FALSE
    ),
    interactions = data.frame(
      from = c("9606.ENSP1", "9606.ENSP1", "9606.ENSP2"),
      to = c("9606.ENSP2", "9606.ENSP3", "9606.ENSP4"),
      combined_score = c(900, 850, 700),
      stringsAsFactors = FALSE
    )
  )

  testthat::expect_true(save_string_network_plot(result, plot_file))
  testthat::expect_true(file.info(plot_file)$size > 1000)
})

testthat::test_that("Immune Deconvolution creates a composition heatmap", {
  testthat::skip_if_not_installed("ggplot2")
  plot_file <- tempfile(fileext = ".png")
  on.exit(unlink(plot_file), add = TRUE)

  result_data <- data.frame(
    cell_type = c("B cell", "CD8 T cell", "Macrophage", "NK cell"),
    Tumor_A = c(0.12, 0.31, 0.40, 0.17),
    Tumor_B = c(0.20, 0.24, 0.33, 0.23),
    Tumor_C = c(0.08, 0.42, 0.28, 0.22),
    check.names = FALSE
  )

  testthat::expect_true(save_immune_composition_plot(result_data, plot_file))
  testthat::expect_true(file.info(plot_file)$size > 1000)
})

testthat::test_that("Immune heatmap separates residual fraction and retains all 50 samples", {
  testthat::skip_if_not_installed("ggplot2")
  sample_names <- paste0("Tumor_", seq_len(50))
  result_data <- data.frame(
    cell_type = c("B cell", "Macrophage M1", "T cell CD8+", "uncharacterized cell"),
    check.names = FALSE
  )
  result_data[sample_names] <- rbind(
    seq(0.01, 0.05, length.out = 50),
    seq(0.03, 0.12, length.out = 50),
    seq(0.02, 0.08, length.out = 50),
    seq(0.90, 0.65, length.out = 50)
  )

  prepared <- prepare_immune_composition_plot(result_data)

  testthat::expect_identical(prepared$total_samples, 50L)
  testthat::expect_identical(prepared$displayed_samples, 50L)
  testthat::expect_false("uncharacterized cell" %in% rownames(prepared$abundance_matrix))
  testthat::expect_equal(nrow(prepared$abundance_matrix), 3L)
  testthat::expect_equal(prepared$residual_summary$median, 0.775)
})

testthat::test_that("Drug Sensitivity creates a ranked response plot", {
  testthat::skip_if_not_installed("ggplot2")
  plot_file <- tempfile(fileext = ".png")
  on.exit(unlink(plot_file), add = TRUE)

  ranking <- data.frame(
    Rank = 1:8,
    Compound = paste("Compound", LETTERS[1:8]),
    Mean_Response = c(0.12, 0.18, 0.24, 0.29, 0.34, 0.41, 0.47, 0.55),
    Measurements = rep(24L, 8)
  )

  testthat::expect_true(
    save_drug_sensitivity_plot(
      ranking,
      plot_file,
      lower_is_sensitive = TRUE
    )
  )
  testthat::expect_true(file.info(plot_file)$size > 1000)
})
