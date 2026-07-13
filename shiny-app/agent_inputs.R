agent_requirements <- list(
  go = list(
    title = "GO Enrichment",
    icon = "GO",
    requirement = "Gene list",
    description = paste(
      "Requires a column containing gene symbols, Entrez IDs,",
      "or another recognized gene identifier."
    ),
    accepted = c(
      "gene", "genes", "gene_symbol", "symbol",
      "entrez", "entrez_id", "gene_id"
    )
  ),

  kegg = list(
    title = "KEGG Pathways",
    icon = "KG",
    requirement = "Gene list",
    description = paste(
      "Requires gene identifiers that can be mapped to",
      "KEGG-supported Entrez gene identifiers."
    ),
    accepted = c(
      "gene", "genes", "gene_symbol", "symbol",
      "entrez", "entrez_id", "gene_id"
    )
  ),

  gsva = list(
    title = "GSVA Activity",
    icon = "GS",
    requirement = "Expression matrix",
    description = paste(
      "Requires genes as rows or columns and multiple",
      "numeric sample-expression columns."
    ),
    accepted = character(0)
  ),

  chea = list(
    title = "ChEA Regulators",
    icon = "TF",
    requirement = "Gene list",
    description = paste(
      "Requires a valid gene-symbol list for",
      "transcription-factor enrichment."
    ),
    accepted = c(
      "gene", "genes", "gene_symbol", "symbol",
      "entrez", "entrez_id", "gene_id"
    )
  )
)

normalize_column_name <- function(value) {
  value |>
    tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("^_|_$", "", x = _)
}

detect_agent_compatibility <- function(data) {
  empty_result <- list(
    dataset_type = "Unknown",
    gene_column = NULL,
    numeric_columns = character(0),
    sample_count = 0L,
    agent_compatibility = c(
      go = FALSE,
      kegg = FALSE,
      gsva = FALSE,
      chea = FALSE
    ),
    messages = c(
      go = "Validate a dataset first.",
      kegg = "Validate a dataset first.",
      gsva = "Validate a dataset first.",
      chea = "Validate a dataset first."
    )
  )

  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    return(empty_result)
  }

  original_names <- names(data)
  normalized_names <- vapply(
    original_names,
    normalize_column_name,
    character(1)
  )

  known_gene_names <- unique(
    unlist(
      lapply(
        agent_requirements[c("go", "kegg", "chea")],
        function(item) item$accepted
      )
    )
  )

  gene_matches <- which(normalized_names %in% known_gene_names)

  gene_column <- if (length(gene_matches) > 0) {
    original_names[[gene_matches[[1]]]]
  } else {
    NULL
  }

  numeric_columns <- original_names[
    vapply(data, is.numeric, logical(1))
  ]

  non_gene_numeric_columns <- setdiff(
    numeric_columns,
    gene_column
  )

  has_gene_list <- !is.null(gene_column)
  has_expression_matrix <- (
    length(non_gene_numeric_columns) >= 2 &&
      nrow(data) >= 5
  )

  dataset_type <- if (has_expression_matrix && has_gene_list) {
    "Gene expression matrix"
  } else if (has_gene_list) {
    "Gene identifier table"
  } else if (has_expression_matrix) {
    "Numeric matrix"
  } else {
    "Unrecognized genomic table"
  }

  compatibility <- c(
    go = has_gene_list,
    kegg = has_gene_list,
    gsva = has_expression_matrix,
    chea = has_gene_list
  )

  gene_message <- if (has_gene_list) {
    paste("Detected gene column:", gene_column)
  } else {
    "No recognized gene identifier column was detected."
  }

  gsva_message <- if (has_expression_matrix) {
    paste(
      length(non_gene_numeric_columns),
      "numeric sample columns detected."
    )
  } else {
    paste(
      "GSVA requires at least two numeric sample columns;",
      length(non_gene_numeric_columns),
      "were detected."
    )
  }

  list(
    dataset_type = dataset_type,
    gene_column = gene_column,
    numeric_columns = non_gene_numeric_columns,
    sample_count = length(non_gene_numeric_columns),
    agent_compatibility = compatibility,
    messages = c(
      go = gene_message,
      kegg = gene_message,
      gsva = gsva_message,
      chea = gene_message
    )
  )
}

agent_selector_card <- function(key) {
  item <- agent_requirements[[key]]

  div(
    id = paste0(key, "-selector-card"),
    class = "agent-selector-card",

    div(
      class = "agent-selector-top",

      div(
        class = paste(
          "agent-selector-icon",
          paste0("selector-icon-", key)
        ),
        item$icon
      ),

      checkboxInput(
        inputId = paste0("enable_", key),
        label = NULL,
        value = TRUE
      )
    ),

    div(
      class = "agent-selector-copy",

      h4(item$title),

      div(
        class = "agent-input-requirement",
        span("INPUT"),
        strong(item$requirement)
      ),

      p(item$description)
    ),

    uiOutput(paste0(key, "_compatibility"))
  )
}

agent_input_control_ui <- function() {
  div(
    class = "agent-input-control glass-card",

    div(
      class = "agent-input-control-heading",

      div(
        div(class = "section-label", "ANALYSIS CONFIGURATION"),
        h2("Choose compatible analysis agents"),
        p(
          class = "section-description",
          paste(
            "The uploaded dataset is inspected automatically.",
            "Only agents compatible with its structure should be selected."
          )
        )
      ),

      uiOutput("selected_agent_summary")
    ),

    uiOutput("detected_dataset_summary"),

    div(
      class = "agent-selector-grid",

      agent_selector_card("go"),
      agent_selector_card("kegg"),
      agent_selector_card("gsva"),
      agent_selector_card("chea")
    )
  )
}
