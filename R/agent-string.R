#' Run STRING protein-protein interaction agent
#'
#' Maps gene symbols to STRING protein identifiers and retrieves interactions
#' between mapped proteins.
#'
#' @param genes Character vector of gene symbols.
#' @param species NCBI taxonomy identifier. Default is 9606 for human.
#' @param score_threshold Minimum STRING combined score from 0 to 1000.
#' @param string_version STRING database version.
#' @param remove_unmapped_rows Logical; remove genes that cannot be mapped.
#'
#' @return A list containing mapped genes, network edges, nodes, hubs, and summary.
#' @export
run_string_agent <- function(
  genes,
  species = 9606,
  score_threshold = 400,
  string_version = "12.0",
  remove_unmapped_rows = TRUE
) {
  if (missing(genes) || length(genes) == 0L) {
    stop("genes must be a non-empty character vector.", call. = FALSE)
  }

  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("STRINGdb is required to run the STRING agent.", call. = FALSE)
  }

  if (
    !is.numeric(score_threshold) ||
    length(score_threshold) != 1L ||
    is.na(score_threshold) ||
    score_threshold < 0 ||
    score_threshold > 1000
  ) {
    stop("score_threshold must be between 0 and 1000.", call. = FALSE)
  }

  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0L) {
    stop("No valid gene symbols remained after preprocessing.", call. = FALSE)
  }

  string_db <- STRINGdb::STRINGdb$new(
    version = string_version,
    species = as.integer(species),
    score_threshold = as.integer(score_threshold),
    input_directory = ""
  )

  gene_table <- data.frame(
    gene_symbol = genes,
    stringsAsFactors = FALSE
  )

  mapped_genes <- suppressMessages(
    string_db$map(
      gene_table,
      "gene_symbol",
      removeUnmappedRows = remove_unmapped_rows,
      takeFirst = TRUE
    )
  )

  if (
    is.null(mapped_genes) ||
    nrow(mapped_genes) == 0L ||
    !"STRING_id" %in% colnames(mapped_genes)
  ) {
    stop("No genes could be mapped to STRING identifiers.", call. = FALSE)
  }

  string_ids <- unique(as.character(mapped_genes$STRING_id))
  interactions <- string_db$get_interactions(string_ids)

  if (is.null(interactions)) {
    interactions <- data.frame()
  }

  if (
    nrow(interactions) > 0L &&
    "combined_score" %in% colnames(interactions)
  ) {
    interactions <- interactions[
      interactions$combined_score >= score_threshold,
      ,
      drop = FALSE
    ]
  }

  nodes <- mapped_genes
  top_hubs <- data.frame()

  if (
    nrow(interactions) > 0L &&
    all(c("from", "to") %in% colnames(interactions))
  ) {
    degree_values <- table(c(interactions$from, interactions$to))

    degree_table <- data.frame(
      STRING_id = names(degree_values),
      degree = as.integer(degree_values),
      stringsAsFactors = FALSE
    )

    nodes <- dplyr::left_join(
      nodes,
      degree_table,
      by = "STRING_id"
    )

    nodes$degree[is.na(nodes$degree)] <- 0L
    nodes <- dplyr::arrange(nodes, dplyr::desc(.data$degree))

    top_hubs <- utils::head(
      nodes[nodes$degree > 0L, , drop = FALSE],
      10L
    )
  } else {
    nodes$degree <- 0L
  }

  summary <- paste0(
    "STRING Agent mapped ",
    nrow(mapped_genes),
    " of ",
    length(genes),
    " input genes and identified ",
    nrow(interactions),
    " protein-protein interactions using a minimum combined score of ",
    score_threshold,
    "."
  )

  list(
    agent_name = "STRING",
    input_genes = genes,
    mapped_genes = mapped_genes,
    nodes = nodes,
    interactions = interactions,
    top_hubs = top_hubs,
    score_threshold = score_threshold,
    summary = summary
  )
}
