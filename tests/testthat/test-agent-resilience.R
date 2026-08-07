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
