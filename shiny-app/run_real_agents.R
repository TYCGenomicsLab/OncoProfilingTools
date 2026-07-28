# Real execution backend for the OncoProfiling Shiny application.
#
# Supported inputs:
# 1. Gene-list / DEG table:
#    GO, KEGG and ChEA
#
# 2. Numeric expression matrix:
#    GSVA
#
# The functions below save results into the paths already used by
# shiny-app/results_helpers.R.

find_project_root <- function() {
  current_directory <- normalizePath(
    getwd(),
    mustWork = TRUE
  )

  if (
    basename(current_directory) == "shiny-app"
  ) {
    return(
      normalizePath(
        file.path(current_directory, ".."),
        mustWork = TRUE
      )
    )
  }

  if (
    dir.exists(
      file.path(current_directory, "shiny-app")
    )
  ) {
    return(current_directory)
  }

  stop(
    paste(
      "Could not locate the OncoProfilingTools project root from:",
      current_directory
    ),
    call. = FALSE
  )
}

real_agent_project_root <- find_project_root()

source(
  file.path(
    real_agent_project_root,
    "R",
    "agent-report.R"
  ),
  local = TRUE
)

source(
  file.path(
    real_agent_project_root,
    "R",
    "agent-go.R"
  ),
  local = TRUE
)

source(
  file.path(
    real_agent_project_root,
    "R",
    "agent-kegg.R"
  ),
  local = TRUE
)

source(
  file.path(
    real_agent_project_root,
    "R",
    "agent-gsva.R"
  ),
  local = TRUE
)

source(
  file.path(
    real_agent_project_root,
    "R",
    "agent-chea.R"
  ),
  local = TRUE
)

source(file.path(real_agent_project_root, "R", "agent-reactome.R"), local = TRUE)
source(file.path(real_agent_project_root, "R", "agent-string.R"), local = TRUE)
source(file.path(real_agent_project_root, "R", "agent-immune.R"), local = TRUE)


# ------------------------------------------------------------
# Input detection helpers
# ------------------------------------------------------------


safe_regex_pattern <- function(pattern, fallback = "") {
  if (
    is.null(pattern) ||
    length(pattern) == 0 ||
    all(is.na(pattern))
  ) {
    return(fallback)
  }

  pattern <- as.character(pattern)
  pattern <- pattern[
    !is.na(pattern) &
      nzchar(pattern)
  ]

  if (length(pattern) == 0) {
    return(fallback)
  }

  pattern[[1]]
}

safe_first_column <- function(candidates, available_columns) {
  candidates <- as.character(candidates)
  available_columns <- as.character(available_columns)

  candidates <- candidates[
    !is.na(candidates) &
      nzchar(candidates)
  ]

  available_columns <- available_columns[
    !is.na(available_columns) &
      nzchar(available_columns)
  ]

  matched <- candidates[
    candidates %in% available_columns
  ]

  if (length(matched) == 0) {
    return(NA_character_)
  }

  matched[[1]]
}

normalise_column_name <- function(value) {
  if (
    is.null(value) ||
    length(value) == 0 ||
    all(is.na(value))
  ) {
    return(character())
  }

  value <- tolower(as.character(value))

  gsub(
    "[^a-z0-9]",
    "",
    value
  )
}


detect_gene_column <- function(data) {
  if (is.null(data) || ncol(data) == 0) {
    return(NULL)
  }

  original_names <- names(data)
  normalised_names <- normalise_column_name(original_names)

  preferred_names <- c(
    "genesymbol",
    "hugosymbol",
    "symbol",
    "gene",
    "genename",
    "geneid",
    "hgncsymbol",
    "externalgenename"
  )

  matched_position <- match(
    preferred_names,
    normalised_names,
    nomatch = 0
  )

  matched_position <- matched_position[
    matched_position > 0
  ]

  if (length(matched_position) > 0) {
    return(
      original_names[
        matched_position[[1]]
      ]
    )
  }

  character_columns <- original_names[
    vapply(
      data,
      function(column) {
        is.character(column) ||
          is.factor(column)
      },
      logical(1)
    )
  ]

  if (length(character_columns) > 0) {
    return(character_columns[[1]])
  }

  NULL
}


clean_gene_symbols <- function(values) {
  genes <- as.character(values)

  genes <- trimws(genes)

  genes <- sub(
    "\\s*\\([^)]*\\)\\s*$",
    "",
    genes
  )

  genes <- sub(
    "\\|.*$",
    "",
    genes
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  genes <- genes[
    !tolower(genes) %in% c(
      "na",
      "nan",
      "null",
      "gene",
      "gene_symbol",
      "symbol"
    )
  ]

  unique(genes)
}


extract_gene_list <- function(data) {
  gene_column <- detect_gene_column(data)

  if (is.null(gene_column)) {
    stop(
      paste(
        "A gene-symbol column could not be detected.",
        "Use a column such as gene_symbol, Gene, SYMBOL or gene."
      ),
      call. = FALSE
    )
  }

  genes <- clean_gene_symbols(
    data[[gene_column]]
  )

  if (length(genes) == 0) {
    stop(
      paste(
        "The detected gene column",
        shQuote(gene_column),
        "does not contain usable gene symbols."
      ),
      call. = FALSE
    )
  }

  list(
    genes = genes,
    column = gene_column
  )
}


detect_expression_matrix <- function(data) {
  if (is.null(data)) {
    return(FALSE)
  }

  if (nrow(data) < 2 || ncol(data) < 10) {
    return(FALSE)
  }

  candidate_data <- data[, -1, drop = FALSE]

  numeric_columns <- vapply(
    candidate_data,
    function(column) {
      converted <- suppressWarnings(
        as.numeric(column)
      )

      valid_fraction <- mean(
        !is.na(converted)
      )

      valid_fraction >= 0.8
    },
    logical(1)
  )

  mean(numeric_columns) >= 0.8
}


prepare_expression_matrix <- function(data) {
  if (!detect_expression_matrix(data)) {
    stop(
      paste(
        "GSVA requires an expression matrix.",
        "Use the first column for sample or gene identifiers",
        "and the remaining columns for numeric expression values."
      ),
      call. = FALSE
    )
  }

  identifiers <- as.character(
    data[[1]]
  )

  expression_data <- data[
    ,
    -1,
    drop = FALSE
  ]

  expression_data[] <- lapply(
    expression_data,
    function(column) {
      suppressWarnings(
        as.numeric(column)
      )
    }
  )

  expression_matrix <- as.matrix(
    expression_data
  )

  storage.mode(expression_matrix) <- "numeric"

  column_gene_fraction <- mean(
    grepl(
      "^[A-Za-z0-9.-]+(?:\\s*\\([^)]*\\))?$",
      colnames(expression_matrix)
    )
  )

  row_gene_fraction <- mean(
    grepl(
      "^[A-Za-z0-9.-]+(?:\\s*\\([^)]*\\))?$",
      identifiers
    )
  )

  if (
    ncol(expression_matrix) > nrow(expression_matrix) ||
      column_gene_fraction >= row_gene_fraction
  ) {
    gene_symbols <- sub(
      "\\s*\\([^)]*\\)\\s*$",
      "",
      colnames(expression_matrix)
    )

    colnames(expression_matrix) <- gene_symbols
    rownames(expression_matrix) <- make.unique(
      identifiers
    )

    expression_matrix <- t(
      expression_matrix
    )
  } else {
    gene_symbols <- sub(
      "\\s*\\([^)]*\\)\\s*$",
      "",
      identifiers
    )

    rownames(expression_matrix) <- gene_symbols
  }

  valid_rows <- (
    !is.na(rownames(expression_matrix)) &
      nzchar(rownames(expression_matrix))
  )

  expression_matrix <- expression_matrix[
    valid_rows,
    ,
    drop = FALSE
  ]

  expression_matrix <- expression_matrix[
    !duplicated(rownames(expression_matrix)),
    ,
    drop = FALSE
  ]

  expression_matrix <- expression_matrix[
    rowSums(is.finite(expression_matrix)) > 0,
    ,
    drop = FALSE
  ]

  expression_matrix
}


# ------------------------------------------------------------
# Plotting helpers
# ------------------------------------------------------------

save_enrichment_dotplot <- function(
  result_data,
  output_file,
  title,
  term_column = "Description"
) {
  if (
    is.null(result_data) ||
      nrow(result_data) == 0 ||
      !term_column %in% names(result_data)
  ) {
    if (file.exists(output_file)) {
      unlink(output_file)
    }

    return(FALSE)
  }

  score_column <- if (
    "p.adjust" %in% names(result_data)
  ) {
    "p.adjust"
  } else if (
    "Adjusted.P.value" %in% names(result_data)
  ) {
    "Adjusted.P.value"
  } else {
    NULL
  }

  if (is.null(score_column)) {
    return(FALSE)
  }

  plot_data <- result_data

  plot_data[[score_column]] <- suppressWarnings(
    as.numeric(plot_data[[score_column]])
  )

  plot_data <- plot_data[
    !is.na(plot_data[[score_column]]) &
      is.finite(plot_data[[score_column]]),
    ,
    drop = FALSE
  ]

  plot_data <- plot_data[
    order(plot_data[[score_column]]),
    ,
    drop = FALSE
  ]

  plot_data <- head(
    plot_data,
    20
  )

  if (nrow(plot_data) == 0) {
    return(FALSE)
  }

  plot_data$display_term <- factor(
    plot_data[[term_column]],
    levels = rev(
      plot_data[[term_column]]
    )
  )

  plot_data$significance <- -log10(
    pmax(
      plot_data[[score_column]],
      .Machine$double.xmin
    )
  )

  if ("Count" %in% names(plot_data)) {
    plot_data$display_count <- suppressWarnings(
      as.numeric(plot_data$Count)
    )
  } else {
    plot_data$display_count <- 3
  }

  plot_data$display_count[
    is.na(plot_data$display_count)
  ] <- 3

  graph <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = significance,
      y = display_term,
      size = display_count,
      color = significance
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.9
    ) +
    ggplot2::labs(
      title = title,
      x = "-log10 adjusted p-value",
      y = NULL,
      size = "Gene count",
      color = "Significance"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 8
      )
    )

  ggplot2::ggsave(
    filename = output_file,
    plot = graph,
    width = 10,
    height = 7,
    dpi = 300
  )

  file.exists(output_file)
}


save_chea_dotplot <- function(
  result_data,
  output_file
) {
  if (
    is.null(result_data) ||
      nrow(result_data) == 0
  ) {
    if (file.exists(output_file)) {
      unlink(output_file)
    }

    return(FALSE)
  }

  required_columns <- c(
    "Term",
    "Adjusted.P.value",
    "Combined.Score"
  )

  if (
    !all(
      required_columns %in% names(result_data)
    )
  ) {
    return(FALSE)
  }

  plot_data <- result_data

  plot_data$Adjusted.P.value <- suppressWarnings(
    as.numeric(
      plot_data$Adjusted.P.value
    )
  )

  plot_data$Combined.Score <- suppressWarnings(
    as.numeric(
      plot_data$Combined.Score
    )
  )

  plot_data <- plot_data[
    !is.na(plot_data$Adjusted.P.value) &
      !is.na(plot_data$Combined.Score),
    ,
    drop = FALSE
  ]

  plot_data <- plot_data[
    order(plot_data$Adjusted.P.value),
    ,
    drop = FALSE
  ]

  plot_data <- head(
    plot_data,
    20
  )

  if (nrow(plot_data) == 0) {
    return(FALSE)
  }

  plot_data$Term <- factor(
    plot_data$Term,
    levels = rev(plot_data$Term)
  )

  plot_data$significance <- -log10(
    pmax(
      plot_data$Adjusted.P.value,
      .Machine$double.xmin
    )
  )

  graph <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Combined.Score,
      y = Term,
      size = significance,
      color = Adjusted.P.value
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.9
    ) +
    ggplot2::labs(
      title = "Top ChEA Transcription-Factor Enrichment",
      x = "Combined score",
      y = NULL,
      size = "-log10 adjusted p-value",
      color = "Adjusted p-value"
    ) +
    ggplot2::theme_classic(
      base_size = 11
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 8
      )
    )

  ggplot2::ggsave(
    filename = output_file,
    plot = graph,
    width = 10,
    height = 7,
    dpi = 300
  )

  file.exists(output_file)
}


# ------------------------------------------------------------
# Individual real agents
# ------------------------------------------------------------

execute_go_agent <- function(
  genes,
  pvalue_cutoff = 0.05
) {
  output_directory <- file.path(
    real_agent_project_root,
    "output",
    "cms4"
  )

  visualization_directory <- file.path(
    real_agent_project_root,
    "output",
    "visualizations"
  )

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    visualization_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  result <- run_go_agent(
    genes = genes,
    ontology = "BP",
    pvalue_cutoff = pvalue_cutoff
  )

  csv_file <- file.path(
    output_directory,
    "go_results.csv"
  )

  plot_file <- file.path(
    visualization_directory,
    "go_biological_process_dotplot.png"
  )

  readr::write_csv(
    result$results,
    csv_file
  )

  save_enrichment_dotplot(
    result_data = result$results,
    output_file = plot_file,
    title = "GO Biological Process Enrichment"
  )

  list(
    success = TRUE,
    agent = "go",
    result = result,
    csv = csv_file,
    plot = plot_file,
    rows = nrow(result$results),
    message = result$summary
  )
}


execute_kegg_agent <- function(
  genes,
  pvalue_cutoff = 0.05
) {
  output_directory <- file.path(
    real_agent_project_root,
    "output",
    "cms4"
  )

  visualization_directory <- file.path(
    real_agent_project_root,
    "output",
    "visualizations"
  )

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    visualization_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  result <- run_kegg_agent(
    genes = genes,
    organism = "hsa",
    pvalue_cutoff = pvalue_cutoff
  )

  csv_file <- file.path(
    output_directory,
    "kegg_results.csv"
  )

  plot_file <- file.path(
    visualization_directory,
    "kegg_pathway_dotplot.png"
  )

  readr::write_csv(
    result$results,
    csv_file
  )

  save_enrichment_dotplot(
    result_data = result$results,
    output_file = plot_file,
    title = "KEGG Pathway Enrichment"
  )

  list(
    success = TRUE,
    agent = "kegg",
    result = result,
    csv = csv_file,
    plot = plot_file,
    rows = nrow(result$results),
    message = result$summary
  )
}

execute_reactome_agent <- function(genes, pvalue_cutoff = 0.05) {
  output_directory <- file.path(real_agent_project_root, "output", "reactome")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  result <- run_reactome_agent(
    genes = genes,
    pvalue_cutoff = pvalue_cutoff
  )
  csv_file <- file.path(output_directory, "reactome_results.csv")
  plot_file <- file.path(output_directory, "reactome_pathways.png")
  readr::write_csv(result$results, csv_file)
  save_enrichment_dotplot(result$results, plot_file, "Reactome Pathway Enrichment")

  list(
    success = TRUE, agent = "reactome", result = result,
    csv = csv_file, plot = plot_file, rows = nrow(result$results),
    message = result$summary
  )
}

execute_wikipathways_agent <- function(genes, pvalue_cutoff = 0.05) {
  stopifnot(requireNamespace("clusterProfiler", quietly = TRUE))
  output_directory <- file.path(real_agent_project_root, "output", "wikipathways")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  result <- clusterProfiler::enrichWP(gene = unique(genes), organism = "Homo sapiens", pvalueCutoff = pvalue_cutoff, qvalueCutoff = 0.2)
  results <- as.data.frame(result)
  csv_file <- file.path(output_directory, "wikipathways_results.csv")
  plot_file <- file.path(output_directory, "wikipathways_pathways.png")
  readr::write_csv(results, csv_file)
  save_enrichment_dotplot(results, plot_file, "WikiPathways Enrichment")
  list(success = TRUE, agent = "wikipathways", result = list(results = results), csv = csv_file, plot = plot_file, rows = nrow(results), message = paste("WikiPathways identified", nrow(results), "enriched pathways."))
}

execute_hallmark_agent <- function(genes, pvalue_cutoff = 0.05) {
  stopifnot(requireNamespace("msigdbr", quietly = TRUE), requireNamespace("clusterProfiler", quietly = TRUE))
  output_directory <- file.path(real_agent_project_root, "output", "hallmark")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  term2gene <- sets[, c("gs_name", "gene_symbol")]
  result <- clusterProfiler::enricher(gene = unique(genes), TERM2GENE = term2gene, pvalueCutoff = pvalue_cutoff, qvalueCutoff = 0.2)
  results <- as.data.frame(result)
  csv_file <- file.path(output_directory, "hallmark_results.csv")
  plot_file <- file.path(output_directory, "hallmark_pathways.png")
  readr::write_csv(results, csv_file)
  save_enrichment_dotplot(results, plot_file, "Cancer Hallmark Enrichment")
  list(success = TRUE, agent = "hallmark", result = list(results = results), csv = csv_file, plot = plot_file, rows = nrow(results), message = paste("Hallmark analysis identified", nrow(results), "enriched gene sets."))
}

execute_string_agent <- function(genes) {
  output_directory <- file.path(real_agent_project_root, "output", "string")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  result <- run_string_agent(genes = genes)
  csv_file <- file.path(output_directory, "string_hub_proteins.csv")
  interaction_file <- file.path(output_directory, "string_interactions.csv")
  readr::write_csv(result$top_hubs, csv_file)
  readr::write_csv(result$interactions, interaction_file)

  list(
    success = TRUE, agent = "string", result = result,
    csv = csv_file, plot = NA_character_, rows = nrow(result$top_hubs),
    message = result$summary
  )
}

execute_immune_agent <- function(data) {
  output_directory <- file.path(real_agent_project_root, "output", "immune")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  expression_matrix <- prepare_expression_matrix(data)
  result <- run_immune_agent(expression_matrix = expression_matrix)
  result_table <- result$results
  csv_file <- file.path(output_directory, "immune_cell_composition.csv")
  readr::write_csv(result_table, csv_file)

  list(
    success = TRUE, agent = "immune", result = result,
    csv = csv_file, plot = NA_character_, rows = nrow(result_table),
    message = result$summary
  )
}

execute_drug_agent <- function(data) {
  original_names <- names(data)
  normalized <- normalise_column_name(original_names)
  if (any(grepl("BRD:", original_names, fixed = TRUE))) {
    numeric_columns <- vapply(data, function(x) mean(!is.na(suppressWarnings(as.numeric(as.character(x))))) >= 0.8, logical(1))
    numeric_columns[[1]] <- FALSE
    values <- data[, numeric_columns, drop = FALSE]
    if (!ncol(values)) stop("No numeric PRISM compound measurements were found.", call. = FALSE)
    ranking <- data.frame(Compound = names(values), Mean_Response = vapply(values, function(x) mean(as.numeric(x), na.rm = TRUE), numeric(1)), Measurements = vapply(values, function(x) sum(is.finite(as.numeric(x))), integer(1)), stringsAsFactors = FALSE)
    ranking <- ranking[order(ranking$Mean_Response, na.last = TRUE), , drop = FALSE]
    ranking$Rank <- seq_len(nrow(ranking))
    ranking <- ranking[, c("Rank", "Compound", "Mean_Response", "Measurements")]
    output_directory <- file.path(real_agent_project_root, "output", "drug")
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    csv_file <- file.path(output_directory, "drug_sensitivity_results.csv")
    readr::write_csv(ranking, csv_file)
    return(list(success = TRUE, agent = "drug", result = list(results = ranking), csv = csv_file, plot = NA_character_, rows = nrow(ranking), message = paste("Drug Sensitivity ranked", nrow(ranking), "PRISM compounds.")))
  }
  compound_candidates <- c("compound", "compoundname", "drug", "drugname", "treatment")
  response_candidates <- c("ic50", "auc", "viability", "sensitivity", "response", "lnic50")
  compound_index <- match(compound_candidates, normalized, nomatch = 0L)
  response_index <- match(response_candidates, normalized, nomatch = 0L)
  compound_index <- compound_index[compound_index > 0L]
  response_index <- response_index[response_index > 0L]

  if (length(compound_index) == 0L || length(response_index) == 0L) {
    stop(
      paste(
        "Drug sensitivity requires a compound/drug column and a numeric",
        "response column such as IC50, AUC, viability, sensitivity, or response."
      ),
      call. = FALSE
    )
  }

  compound_column <- original_names[compound_index[[1]]]
  response_column <- original_names[response_index[[1]]]
  response_values <- suppressWarnings(as.numeric(data[[response_column]]))
  valid <- !is.na(data[[compound_column]]) & nzchar(trimws(as.character(data[[compound_column]]))) & is.finite(response_values)
  response_data <- data.frame(
    Compound = trimws(as.character(data[[compound_column]][valid])),
    Response = response_values[valid],
    stringsAsFactors = FALSE
  )
  if (nrow(response_data) == 0L) stop("No valid drug-response measurements were found.", call. = FALSE)

  summary_table <- stats::aggregate(Response ~ Compound, response_data, function(x) c(mean = mean(x), n = length(x)))
  ranking <- data.frame(
    Compound = summary_table$Compound,
    Mean_Response = summary_table$Response[, "mean"],
    Measurements = as.integer(summary_table$Response[, "n"]),
    stringsAsFactors = FALSE
  )
  lower_is_sensitive <- grepl("ic50|auc|viability", normalise_column_name(response_column))
  ranking <- ranking[order(ranking$Mean_Response, decreasing = !lower_is_sensitive), , drop = FALSE]
  ranking$Rank <- seq_len(nrow(ranking))
  ranking <- ranking[, c("Rank", "Compound", "Mean_Response", "Measurements")]

  output_directory <- file.path(real_agent_project_root, "output", "drug")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  csv_file <- file.path(output_directory, "drug_sensitivity_results.csv")
  readr::write_csv(ranking, csv_file)

  list(
    success = TRUE, agent = "drug",
    result = list(
      agent_name = "Drug Sensitivity", results = ranking,
      response_column = response_column, lower_is_sensitive = lower_is_sensitive
    ),
    csv = csv_file, plot = NA_character_, rows = nrow(ranking),
    message = paste("Drug Sensitivity ranked", nrow(ranking), "compounds using", response_column, "measurements.")
  )
}


execute_chea_agent <- function(
  genes
) {
  output_directory <- file.path(
    real_agent_project_root,
    "output",
    "chea_cms4"
  )

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Attach enrichR so its connection options are initialized.
  suppressPackageStartupMessages(
    library(enrichR)
  )

  # Human Enrichr is the default site after package initialization.
  # Avoid setEnrichrSite() here because background workers may not
  # have enrichR.sites.base.address initialized before attachment.
  result <- run_chea_agent(
    genes = genes,
    database = "ChEA_2022"
  )

  csv_file <- file.path(
    output_directory,
    "chea_results.csv"
  )

  plot_file <- file.path(
    output_directory,
    "chea_tf_dotplot.png"
  )

  readr::write_csv(
    result$results,
    csv_file
  )

  save_chea_dotplot(
    result_data = result$results,
    output_file = plot_file
  )

  list(
    success = TRUE,
    agent = "chea",
    result = result,
    csv = csv_file,
    plot = plot_file,
    rows = nrow(result$results),
    message = result$summary
  )
}


execute_gsva_agent <- function(data) {
  output_directory <- file.path(
    real_agent_project_root,
    "output",
    "gsva_bowel"
  )

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  expression_matrix <- prepare_expression_matrix(
    data
  )

  hallmark_data <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "H"
  )

  hallmark_sets <- split(
    hallmark_data$gene_symbol,
    hallmark_data$gs_name
  )

  result <- run_gsva_agent(
    expression_matrix = expression_matrix,
    gene_sets = hallmark_sets,
    kcdf = "Gaussian",
    min_size = 10,
    max_size = 500
  )

  csv_file <- file.path(
    output_directory,
    "gsva_hallmark_scores.csv"
  )

  plot_file <- file.path(
    output_directory,
    "gsva_hallmark_heatmap.png"
  )

  result_table <- as.data.frame(
    result$results
  )

  result_table <- tibble::rownames_to_column(
    result_table,
    var = "Pathway"
  )

  readr::write_csv(
    result_table,
    csv_file
  )

  pathway_variance <- apply(
    result$results,
    1,
    stats::var,
    na.rm = TRUE
  )

  pathway_variance <- pathway_variance[
    is.finite(pathway_variance)
  ]

  top_count <- min(
    20,
    length(pathway_variance)
  )

  if (top_count > 0) {
    top_pathways <- names(
      sort(
        pathway_variance,
        decreasing = TRUE
      )
    )[seq_len(top_count)]

    heatmap_matrix <- result$results[
      top_pathways,
      ,
      drop = FALSE
    ]

    grDevices::png(
      filename = plot_file,
      width = 2400,
      height = 1800,
      res = 220
    )

    pheatmap::pheatmap(
      heatmap_matrix,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      main = "Top Variable Hallmark Pathway Activity",
      border_color = NA
    )

    grDevices::dev.off()
  }

  list(
    success = TRUE,
    agent = "gsva",
    result = result,
    csv = csv_file,
    plot = plot_file,
    rows = nrow(result_table),
    message = result$summary
  )
}


# ------------------------------------------------------------
# Safe dispatcher used by the Shiny server
# ------------------------------------------------------------

run_selected_real_agent <- function(
  agent,
  data,
  pvalue_cutoff = 0.05
) {
  agent <- tolower(agent)

  tryCatch(
    {
      if (agent == "gsva") {
        return(
          execute_gsva_agent(data)
        )
      }

      if (agent == "immune") return(execute_immune_agent(data))
      if (agent == "drug") return(execute_drug_agent(data))

      gene_input <- extract_gene_list(
        data
      )

      if (length(gene_input$genes) < 2) {
        stop(
          paste(
            "At least two usable gene symbols are required",
            "for enrichment analysis."
          ),
          call. = FALSE
        )
      }

      result <- switch(
        agent,

        go = execute_go_agent(
          genes = gene_input$genes,
          pvalue_cutoff = pvalue_cutoff
        ),

        kegg = execute_kegg_agent(
          genes = gene_input$genes,
          pvalue_cutoff = pvalue_cutoff
        ),

        reactome = execute_reactome_agent(
          genes = gene_input$genes,
          pvalue_cutoff = pvalue_cutoff
        ),

        wikipathways = execute_wikipathways_agent(genes = gene_input$genes, pvalue_cutoff = pvalue_cutoff),

        string = execute_string_agent(
          genes = gene_input$genes
        ),

        hallmark = execute_hallmark_agent(genes = gene_input$genes, pvalue_cutoff = pvalue_cutoff),

        chea = execute_chea_agent(
          genes = gene_input$genes
        ),

        stop(
          paste(
            "Unknown analysis agent:",
            agent
          ),
          call. = FALSE
        )
      )

      result$gene_column <- gene_input$column
      result$input_gene_count <- length(
        gene_input$genes
      )

      result
    },
    error = function(error) {
      list(
        success = FALSE,
        agent = agent,
        rows = 0L,
        message = conditionMessage(error),
        error = conditionMessage(error)
      )
    }
  )
}


# Final GSVA orientation override:
# Always returns genes in rows and samples in columns.
prepare_expression_matrix <- function(data, ...) {

  data <- as.data.frame(
    data,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("Expression dataset is empty.")
  }

  original_names <- names(data)

  normalized_names <- tolower(
    gsub(
      "[^a-z0-9]+",
      "_",
      original_names
    )
  )

  gene_candidates <- c(
    "gene_symbol",
    "gene",
    "genes",
    "symbol",
    "gene_name",
    "hgnc_symbol"
  )

  gene_index <- which(
    normalized_names %in% gene_candidates
  )

  numeric_index <- which(
    vapply(
      data,
      function(column) {
        is.numeric(column) ||
          sum(
            !is.na(
              suppressWarnings(
                as.numeric(as.character(column))
              )
            )
          ) >= max(5, floor(length(column) * 0.8))
      },
      logical(1)
    )
  )

  if (length(gene_index) > 0) {

    gene_index <- gene_index[[1]]

    sample_index <- setdiff(
      numeric_index,
      gene_index
    )

    if (length(sample_index) < 2) {
      stop(
        "GSVA requires at least two numeric sample columns."
      )
    }

    genes <- as.character(
      data[[gene_index]]
    )

    expression_matrix <- as.matrix(
      data[, sample_index, drop = FALSE]
    )

    storage.mode(expression_matrix) <- "numeric"

    rownames(expression_matrix) <- genes

  } else {

    # Handle a matrix stored as samples × genes.
    expression_matrix <- as.matrix(data)
    storage.mode(expression_matrix) <- "numeric"

    if (nrow(expression_matrix) < ncol(expression_matrix)) {
      expression_matrix <- t(expression_matrix)
    }
  }

  genes <- rownames(expression_matrix)

  if (
    is.null(genes) ||
    all(grepl("^[0-9]+$", genes))
  ) {
    stop(
      paste(
        "GSVA could not find gene symbols.",
        "The gene-symbol column must be named gene_symbol,",
        "gene, symbol, gene_name, or hgnc_symbol."
      )
    )
  }

  genes <- trimws(
    sub(
      "\\s*\\([^)]*\\)\\s*$",
      "",
      genes
    )
  )

  genes <- toupper(genes)

  valid_rows <- (
    !is.na(genes) &
      nzchar(genes) &
      rowSums(is.finite(expression_matrix)) > 0
  )

  expression_matrix <- expression_matrix[
    valid_rows,
    ,
    drop = FALSE
  ]

  genes <- genes[valid_rows]

  # Average duplicated genes instead of adding .1/.2 suffixes.
  if (anyDuplicated(genes)) {

    summed <- rowsum(
      expression_matrix,
      group = genes,
      reorder = FALSE,
      na.rm = TRUE
    )

    counts <- table(genes)[rownames(summed)]

    expression_matrix <- summed /
      as.numeric(counts)

  } else {
    rownames(expression_matrix) <- genes
  }

  expression_matrix[
    !is.finite(expression_matrix)
  ] <- NA_real_

  keep_rows <- rowSums(
    !is.na(expression_matrix)
  ) >= 2

  expression_matrix <- expression_matrix[
    keep_rows,
    ,
    drop = FALSE
  ]

  if (nrow(expression_matrix) < 10) {
    stop(
      "Too few valid gene rows remain for GSVA."
    )
  }

  if (ncol(expression_matrix) < 2) {
    stop(
      "Too few numeric sample columns remain for GSVA."
    )
  }

  expression_matrix
}
