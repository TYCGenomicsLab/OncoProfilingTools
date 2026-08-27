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


looks_like_gene_identifier <- function(values) {
  values <- trimws(as.character(values))
  valid <- !is.na(values) & nzchar(values)
  result <- grepl("^ENSG[0-9]{11}(?:[.][0-9]+)?$", values, ignore.case = TRUE) |
    grepl("^[0-9]{1,12}$", values) |
    grepl("^[A-Za-z][A-Za-z0-9.-]{1,30}$", values)
  result[!valid] <- FALSE
  result
}


detect_gene_identifier_type <- function(values) {
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return("unknown")

  ensembl <- mean(grepl("^ENSG[0-9]{11}(?:[.][0-9]+)?$", values, ignore.case = TRUE))
  entrez <- mean(grepl("^[0-9]{1,12}$", values))
  symbols <- mean(grepl("^[A-Za-z][A-Za-z0-9.-]{1,30}$", values))

  if (ensembl >= 0.8) return("ensembl")
  if (entrez >= 0.8) return("entrez")
  if (symbols >= 0.8) return("symbol")
  if (max(ensembl, entrez, symbols) >= 0.3) return("mixed")
  "unknown"
}


read_analysis_dataset <- function(path, name) {
  extension <- tolower(tools::file_ext(name))
  delimiter <- switch(extension, csv = ",", tsv = "\t", txt = "\t", NULL)
  if (is.null(delimiter)) {
    stop("Upload a CSV, TSV, or TXT file.", call. = FALSE)
  }

  preview_lines <- readLines(path, n = 8L, warn = FALSE, encoding = "UTF-8")
  preview_lines <- preview_lines[nzchar(trimws(preview_lines))]
  if (!length(preview_lines)) stop("The uploaded file is empty.", call. = FALSE)

  if (identical(extension, "txt") && !any(grepl("\t", preview_lines, fixed = TRUE))) {
    delimiter <- if (any(grepl(",", preview_lines, fixed = TRUE))) "," else "\t"
  }

  split_line <- function(line) {
    if (identical(delimiter, "\t") && !grepl("\t", line, fixed = TRUE)) {
      return(trimws(line))
    }
    trimws(strsplit(line, delimiter, fixed = TRUE)[[1L]])
  }

  first_fields <- split_line(preview_lines[[1L]])
  normalized_first <- normalise_column_name(first_fields)
  recognized_headers <- c(
    "genesymbol", "hugosymbol", "hgncsymbol", "gene", "genes", "symbol",
    "geneid", "ensembl", "ensemblgeneid", "entrezid", "genename", "genetype",
    "compound", "compoundname", "drug", "drugname", "treatment", "response",
    "auc", "ic50", "viability", "sensitivity", "log2foldchange", "logfc",
    "pvalue", "padj", "fdr", "qvalue"
  )
  explicit_header <- any(normalized_first %in% recognized_headers)
  data_like_first <- mean(
    looks_like_gene_identifier(first_fields) |
      grepl("^-?[0-9]+(?:[.][0-9]+)?(?:e[+-]?[0-9]+)?$", first_fields, ignore.case = TRUE)
  ) >= 0.8
  header_detected <- explicit_header || !data_like_first

  column_names <- if (header_detected) TRUE else FALSE
  data <- readr::read_delim(
    path,
    delim = delimiter,
    col_names = column_names,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "unique",
    trim_ws = TRUE
  )

  if (!header_detected) {
    names(data) <- if (ncol(data) == 1L) "gene_id" else paste0("column_", seq_len(ncol(data)))
  }

  attr(data, "import_metadata") <- list(
    source_name = basename(name),
    extension = extension,
    delimiter = if (identical(delimiter, "\t")) "tab" else "comma",
    header_detected = header_detected
  )
  data
}


map_gene_identifiers <- function(values, identifier_type = NULL) {
  input <- trimws(as.character(values))
  input <- input[!is.na(input) & nzchar(input)]
  input <- unique(input)
  if (is.null(identifier_type) || !length(identifier_type) || is.na(identifier_type[[1L]])) {
    identifier_type <- detect_gene_identifier_type(input)
  }
  identifier_type <- as.character(identifier_type[[1L]])

  if (!length(input)) {
    return(list(
      symbols = character(),
      entrez_ids = character(),
      metadata = list(
        identifier_type = identifier_type,
        input_count = 0L,
        mapped_count = 0L,
        unmapped_count = 0L,
        duplicate_mappings_removed = 0L,
        mapping_rate = 0,
        unmapped_examples = character()
      ),
      mapping_table = data.frame()
    ))
  }

  if (identical(identifier_type, "symbol")) {
    symbols <- toupper(input)
    table <- data.frame(
      input_id = input,
      normalized_id = symbols,
      SYMBOL = symbols,
      ENTREZID = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop(
        "Human identifier mapping requires AnnotationDbi and org.Hs.eg.db.",
        call. = FALSE
      )
    }

    key_type <- switch(
      identifier_type,
      ensembl = "ENSEMBL",
      entrez = "ENTREZID",
      mixed = "ENSEMBL",
      NULL
    )
    if (is.null(key_type)) {
      stop(
        "Gene identifiers were not recognized. Use HGNC symbols, Ensembl gene IDs, or Entrez IDs.",
        call. = FALSE
      )
    }
    normalized <- if (identical(key_type, "ENSEMBL")) {
      toupper(sub("[.][0-9]+$", "", input))
    } else {
      sub("[.]0+$", "", input)
    }
    selected <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(normalized),
      keytype = key_type,
      columns = c("SYMBOL", "ENTREZID")
    ))
    names(selected)[names(selected) == key_type] <- "normalized_id"
    table <- merge(
      data.frame(input_id = input, normalized_id = normalized, stringsAsFactors = FALSE),
      selected,
      by = "normalized_id",
      all.x = TRUE,
      sort = FALSE
    )
    table <- table[match(input, table$input_id), , drop = FALSE]
  }

  mapped <- !is.na(table$SYMBOL) & nzchar(trimws(as.character(table$SYMBOL)))
  symbols_with_duplicates <- toupper(trimws(as.character(table$SYMBOL[mapped])))
  symbols <- unique(symbols_with_duplicates)
  entrez_ids <- unique(as.character(table$ENTREZID[mapped]))
  entrez_ids <- entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)]
  mapped_input_count <- length(unique(table$input_id[mapped]))

  list(
    symbols = symbols,
    entrez_ids = entrez_ids,
    mapping_table = table,
    metadata = list(
      identifier_type = identifier_type,
      input_count = length(input),
      mapped_count = mapped_input_count,
      output_symbol_count = length(symbols),
      unmapped_count = length(input) - mapped_input_count,
      duplicate_mappings_removed = length(symbols_with_duplicates) - length(symbols),
      mapping_rate = mapped_input_count / length(input),
      unmapped_examples = utils::head(table$input_id[!mapped], 8L)
    )
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
    "ensembl",
    "ensemblgeneid",
    "entrezid",
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


detect_analysis_column <- function(data, candidates, pattern = NULL) {
  normalized <- normalise_column_name(names(data))
  position <- match(candidates, normalized, nomatch = 0L)
  position <- position[position > 0L]
  if (length(position)) return(names(data)[position[[1L]]])
  if (!is.null(pattern)) {
    position <- which(grepl(pattern, normalized))
    if (length(position)) return(names(data)[position[[1L]]])
  }
  NULL
}


prepare_gene_input <- function(
  data,
  pvalue_cutoff = 0.05,
  effect_cutoff = 1,
  max_genes = 2000L,
  max_unranked_genes = 5000L
) {
  gene_column <- detect_gene_column(data)

  if (is.null(gene_column)) {
    stop(
      paste(
        "A gene identifier column could not be detected.",
        "Use HGNC symbols, Ensembl gene IDs, or Entrez IDs in a named or single-column file."
      ),
      call. = FALSE
    )
  }

  raw_genes <- trimws(as.character(data[[gene_column]]))
  usable <- !is.na(raw_genes) & nzchar(raw_genes) &
    !tolower(raw_genes) %in% c("na", "nan", "null", "gene", "gene_symbol", "symbol")
  original_gene_count <- length(unique(raw_genes[usable]))

  if (original_gene_count == 0L) {
    stop(
      paste(
        "The detected gene column",
        shQuote(gene_column),
        "does not contain usable gene identifiers."
      ),
      call. = FALSE
    )
  }

  adjusted_p_column <- detect_analysis_column(
    data,
    c("padj", "adjustedpvalue", "adjustedp", "adjpvalue", "adjp", "fdr", "qvalue"),
    pattern = "^(padj|fdr|qvalue)|adjusted.*pvalue"
  )
  if (is.null(adjusted_p_column)) {
    adjusted_p_column <- detect_analysis_column(
      data,
      c("pvalue", "pval", "pvalue2sided"),
      pattern = "^pvalue$|^pval$"
    )
  }
  effect_column <- detect_analysis_column(
    data,
    c(
      "log2foldchange", "logfoldchange", "log2fc", "logfc",
      "avglog2fc", "avglogfc", "difference", "diff", "effectsize"
    ),
    pattern = "log2?.*fold.*change|(^|mean)diff|diff.*(minus|vs)|effectsize"
  )

  selected <- usable
  p_values <- NULL
  effects <- NULL
  selection_parts <- character()

  if (!is.null(adjusted_p_column)) {
    p_values <- suppressWarnings(as.numeric(as.character(data[[adjusted_p_column]])))
    selected <- selected & is.finite(p_values) & p_values <= pvalue_cutoff
    selection_parts <- c(
      selection_parts,
      paste0(adjusted_p_column, " ≤ ", format(pvalue_cutoff, trim = TRUE))
    )
  }

  if (!is.null(effect_column)) {
    effects <- suppressWarnings(as.numeric(as.character(data[[effect_column]])))
    base_selected <- selected
    selected_with_effect <- base_selected & is.finite(effects) & abs(effects) >= effect_cutoff

    if (sum(selected_with_effect) < 2L && !is.null(p_values)) {
      selected_with_effect <- base_selected & is.finite(effects) & abs(effects) >= 0.5
      if (sum(selected_with_effect) >= 2L) effect_cutoff <- 0.5
    }
    if (sum(selected_with_effect) < 2L && !is.null(p_values)) {
      selected_with_effect <- base_selected
      effect_column <- NULL
      effects <- NULL
    }
    selected <- selected_with_effect
    if (!is.null(effect_column)) {
      selection_parts <- c(
        selection_parts,
        paste0("|", effect_column, "| ≥ ", format(effect_cutoff, trim = TRUE))
      )
    }
  }

  has_statistics <- !is.null(adjusted_p_column) || !is.null(effect_column)
  if (!has_statistics && original_gene_count > max_unranked_genes) {
    stop(
      paste0(
        "The uploaded table contains ", format(original_gene_count, big.mark = ","),
        " genes but no p-value or effect-size column was detected. ",
        "A near-whole-genome list is not a meaningful over-representation input. ",
        "Upload a filtered gene list or include columns such as adjusted p-value and log2 fold change."
      ),
      call. = FALSE
    )
  }

  selected_rows <- which(selected)
  if (!length(selected_rows)) {
    stop(
      paste0(
        "No genes passed the automatic analysis filter",
        if (length(selection_parts)) paste0(" (", paste(selection_parts, collapse = "; "), ")") else "",
        ". Review the input statistics or upload a curated gene list."
      ),
      call. = FALSE
    )
  }

  if (!is.null(p_values)) {
    selected_rows <- selected_rows[order(
      p_values[selected_rows],
      if (!is.null(effects)) -abs(effects[selected_rows]) else seq_along(selected_rows),
      na.last = TRUE
    )]
  } else if (!is.null(effects)) {
    selected_rows <- selected_rows[order(-abs(effects[selected_rows]), na.last = TRUE)]
  }

  selected_rows <- utils::head(selected_rows, max(2L, as.integer(max_genes)))
  selected_genes <- raw_genes[selected_rows]
  unique_rows <- !duplicated(selected_genes)
  selected_rows <- selected_rows[unique_rows]
  selected_genes <- selected_genes[unique_rows]

  if (length(selected_genes) < 2L) {
    stop("At least two genes must pass the analysis filter.", call. = FALSE)
  }

  mapping <- map_gene_identifiers(selected_genes)
  if (length(mapping$symbols) < 2L) {
    stop(
      paste0(
        "Only ", length(mapping$symbols), " identifier",
        if (length(mapping$symbols) == 1L) "" else "s",
        " mapped to a current human HGNC symbol. At least two mapped genes are required."
      ),
      call. = FALSE
    )
  }

  method <- if (length(selection_parts)) {
    paste(selection_parts, collapse = " and ")
  } else {
    "the supplied curated gene list"
  }
  capped <- sum(selected) > length(selected_genes)
  selection_note <- paste0(
    "Analysis selected ", format(length(selected_genes), big.mark = ","), " of ",
    format(original_gene_count, big.mark = ","), " usable genes using ", method,
    if (capped) paste0(" (capped at the ", format(length(selected_genes), big.mark = ","), " strongest rows)") else "",
    if (!is.null(effects)) "; positive and negative effects are combined" else "",
    ". Identifier type: ", mapping$metadata$identifier_type,
    "; mapped ", mapping$metadata$mapped_count, "/", mapping$metadata$input_count,
    " inputs to ", mapping$metadata$output_symbol_count, " unique HGNC symbols."
  )

  list(
    genes = mapping$symbols,
    entrez_ids = mapping$entrez_ids,
    selected_identifiers = selected_genes,
    mapping = mapping$mapping_table,
    mapping_metadata = mapping$metadata,
    column = gene_column,
    data = data[selected_rows, , drop = FALSE],
    original_gene_count = original_gene_count,
    selection_note = selection_note,
    pvalue_column = adjusted_p_column,
    effect_column = effect_column
  )
}


extract_gene_list <- function(data) {
  prepare_gene_input(data)
}


clean_wikipathways_entrez_ids <- function(values) {
  values <- trimws(as.character(values))
  values <- sub("\\.0+$", "", values)
  values <- values[
    !is.na(values) &
      nzchar(values) &
      grepl("^[0-9]+$", values) &
      values != "0"
  ]
  unique(values)
}


extract_wikipathways_entrez_ids <- function(data) {
  if (is.null(data) || ncol(data) == 0L) {
    return(character())
  }

  normalized_names <- normalise_column_name(names(data))
  candidates <- c(
    "entrezid",
    "entrezgeneid",
    "ncbigeneid",
    "ncbigene"
  )
  position <- match(candidates, normalized_names, nomatch = 0L)
  position <- position[position > 0L]

  if (length(position) == 0L) {
    return(character())
  }

  clean_wikipathways_entrez_ids(data[[position[[1L]]]])
}


map_wikipathways_symbols_to_entrez <- function(genes) {
  if (
    !requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("RSQLite", quietly = TRUE)
  ) {
    stop(
      "WikiPathways symbol mapping requires the DBI and RSQLite packages.",
      call. = FALSE
    )
  }

  database_file <- system.file(
    "extdata",
    "org.Hs.eg.sqlite",
    package = "org.Hs.eg.db"
  )

  if (!nzchar(database_file) || !file.exists(database_file)) {
    stop(
      "The org.Hs.eg.db annotation database is unavailable.",
      call. = FALSE
    )
  }

  connection <- DBI::dbConnect(
    RSQLite::SQLite(),
    database_file
  )
  on.exit(DBI::dbDisconnect(connection), add = TRUE)

  chunks <- split(
    genes,
    ceiling(seq_along(genes) / 900L)
  )
  mappings <- lapply(chunks, function(chunk) {
    quoted_symbols <- DBI::dbQuoteString(connection, chunk)
    statement <- paste0(
      "SELECT gi.symbol, g.gene_id ",
      "FROM gene_info AS gi ",
      "JOIN genes AS g USING (_id) ",
      "WHERE gi.symbol IN (",
      paste(quoted_symbols, collapse = ","),
      ")"
    )
    DBI::dbGetQuery(connection, statement)
  })

  mapped <- do.call(rbind, mappings)
  if (is.null(mapped) || nrow(mapped) == 0L) {
    return(character())
  }

  clean_wikipathways_entrez_ids(mapped$gene_id)
}


fetch_wikipathways_gmt <- function(cache_directory) {
  base_url <- "https://data.wikipathways.org/current/gmt/"
  pattern <- "wikipathways-[0-9]+-gmt-Homo_sapiens\\.gmt"

  dir.create(
    cache_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  cached_files <- list.files(
    cache_directory,
    pattern = paste0("^", pattern, "$"),
    full.names = TRUE
  )
  cached_files <- cached_files[
    file.info(cached_files)$size > 0
  ]
  cached_files <- cached_files[
    order(basename(cached_files), decreasing = TRUE)
  ]

  if (
    length(cached_files) > 0L &&
      difftime(
        Sys.time(),
        file.info(cached_files[[1L]])$mtime,
        units = "days"
      ) < 7
  ) {
    return(cached_files[[1L]])
  }

  index <- tryCatch(
    readLines(base_url, warn = FALSE),
    error = function(error) {
      NULL
    }
  )

  if (is.null(index)) {
    if (length(cached_files) > 0L) {
      return(cached_files[[1L]])
    }
    stop(
      "Could not read the WikiPathways GMT index and no cache is available.",
      call. = FALSE
    )
  }

  matches <- regmatches(
    index,
    gregexpr(pattern, index, perl = TRUE)
  )
  filenames <- sort(
    unique(unlist(matches, use.names = FALSE)),
    decreasing = TRUE
  )

  if (length(filenames) == 0L) {
    if (length(cached_files) > 0L) {
      return(cached_files[[1L]])
    }
    stop(
      "The current human WikiPathways GMT file could not be located.",
      call. = FALSE
    )
  }

  cache_file <- file.path(
    cache_directory,
    filenames[[1L]]
  )

  if (!file.exists(cache_file) || file.info(cache_file)$size == 0) {
    temporary_file <- tempfile(
      "wikipathways-",
      tmpdir = cache_directory,
      fileext = ".gmt"
    )
    on.exit(unlink(temporary_file), add = TRUE)

    tryCatch(
      utils::download.file(
        paste0(base_url, filenames[[1L]]),
        temporary_file,
        mode = "wb",
        quiet = TRUE
      ),
      error = function(error) {
        stop(
          paste(
            "Could not download the WikiPathways GMT file:",
            conditionMessage(error)
          ),
          call. = FALSE
        )
      }
    )

    if (
      !file.rename(temporary_file, cache_file) &&
        !file.exists(cache_file)
    ) {
      stop(
        "The WikiPathways GMT file could not be cached.",
        call. = FALSE
      )
    }
  }

  cache_file
}


calculate_wikipathways_enrichment <- function(
  entrez_ids,
  gmt_file,
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.20,
  min_gene_set_size = 10L,
  max_gene_set_size = 500L
) {
  if (!requireNamespace("qvalue", quietly = TRUE)) {
    stop(
      "WikiPathways enrichment requires the qvalue package.",
      call. = FALSE
    )
  }

  lines <- readLines(gmt_file, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 3L]

  metadata <- strsplit(
    vapply(fields, `[[`, character(1), 1L),
    "%",
    fixed = TRUE
  )
  valid <- lengths(metadata) >= 3L
  fields <- fields[valid]
  metadata <- metadata[valid]

  pathway_ids <- vapply(
    metadata,
    `[[`,
    character(1),
    3L
  )
  pathway_names <- vapply(
    metadata,
    `[[`,
    character(1),
    1L
  )
  gene_sets <- lapply(fields, function(values) {
    unique(values[-c(1L, 2L)])
  })

  universe <- unique(unlist(gene_sets, use.names = FALSE))
  input_genes <- intersect(
    clean_wikipathways_entrez_ids(entrez_ids),
    universe
  )

  if (length(input_genes) == 0L) {
    stop(
      paste(
        "None of the supplied genes mapped to the current human",
        "WikiPathways collection."
      ),
      call. = FALSE
    )
  }

  gene_set_sizes <- lengths(gene_sets)
  tested <- gene_set_sizes >= min_gene_set_size &
    gene_set_sizes <= max_gene_set_size
  pathway_ids <- pathway_ids[tested]
  pathway_names <- pathway_names[tested]
  gene_sets <- gene_sets[tested]
  gene_set_sizes <- gene_set_sizes[tested]

  overlaps <- lapply(gene_sets, intersect, y = input_genes)
  counts <- lengths(overlaps)
  pvalues <- stats::phyper(
    counts - 1L,
    gene_set_sizes,
    length(universe) - gene_set_sizes,
    length(input_genes),
    lower.tail = FALSE
  )
  adjusted_pvalues <- stats::p.adjust(
    pvalues,
    method = "BH"
  )
  qvalues <- tryCatch(
    qvalue::qvalue(pvalues)$qvalues,
    error = function(error) {
      rep(NA_real_, length(pvalues))
    }
  )

  results <- data.frame(
    ID = pathway_ids,
    Description = pathway_names,
    GeneRatio = paste0(counts, "/", length(input_genes)),
    BgRatio = paste0(gene_set_sizes, "/", length(universe)),
    pvalue = pvalues,
    p.adjust = adjusted_pvalues,
    qvalue = qvalues,
    geneID = vapply(
      overlaps,
      paste,
      collapse = "/",
      FUN.VALUE = character(1)
    ),
    Count = counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  significant <- results$pvalue <= pvalue_cutoff &
    results$p.adjust <= pvalue_cutoff &
    (is.na(results$qvalue) | results$qvalue <= qvalue_cutoff)
  results <- results[significant, , drop = FALSE]
  results <- results[order(results$pvalue), , drop = FALSE]
  rownames(results) <- NULL

  list(
    results = results,
    mapped_genes = input_genes
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

  display_labels <- gsub("_", " ", as.character(plot_data[[term_column]]), fixed = TRUE)
  display_labels <- sub("^HALLMARK ", "", display_labels)
  display_labels <- vapply(
    display_labels,
    function(label) paste(strwrap(label, width = 52), collapse = "\n"),
    character(1)
  )
  plot_data$display_term <- factor(
    display_labels,
    levels = rev(unique(display_labels))
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


save_string_network_plot <- function(result, output_file, max_nodes = 35L) {
  nodes <- as.data.frame(if (is.null(result$nodes)) data.frame() else result$nodes, stringsAsFactors = FALSE)
  interactions <- as.data.frame(if (is.null(result$interactions)) data.frame() else result$interactions, stringsAsFactors = FALSE)

  required_node_columns <- c("STRING_id", "gene_symbol", "degree")
  if (!nrow(nodes) || !all(required_node_columns %in% names(nodes))) {
    if (file.exists(output_file)) unlink(output_file)
    return(FALSE)
  }

  nodes$degree <- suppressWarnings(as.numeric(nodes$degree))
  nodes <- nodes[is.finite(nodes$degree) & nodes$degree > 0, , drop = FALSE]
  nodes <- nodes[order(nodes$degree, decreasing = TRUE), , drop = FALSE]
  nodes <- utils::head(nodes, as.integer(max_nodes))
  if (!nrow(nodes)) {
    if (file.exists(output_file)) unlink(output_file)
    return(FALSE)
  }

  can_draw_network <- requireNamespace("igraph", quietly = TRUE) &&
    nrow(interactions) > 0L && all(c("from", "to") %in% names(interactions))

  if (can_draw_network) {
    selected_ids <- as.character(nodes$STRING_id)
    edges <- interactions[
      as.character(interactions$from) %in% selected_ids &
        as.character(interactions$to) %in% selected_ids,
      ,
      drop = FALSE
    ]

    if (nrow(edges)) {
      edge_table <- data.frame(
        from = as.character(edges$from),
        to = as.character(edges$to),
        stringsAsFactors = FALSE
      )
      vertices <- data.frame(
        name = as.character(nodes$STRING_id),
        label = as.character(nodes$gene_symbol),
        degree = nodes$degree,
        stringsAsFactors = FALSE
      )
      graph <- igraph::graph_from_data_frame(edge_table, directed = FALSE, vertices = vertices)
      graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

      grDevices::png(output_file, width = 2400, height = 1700, res = 220, bg = "white")
      set.seed(42L)
      layout <- igraph::layout_with_fr(graph)
      degree_values <- igraph::vertex_attr(graph, "degree")
      label_values <- igraph::vertex_attr(graph, "label")
      size_values <- 8 + 11 * sqrt(degree_values / max(degree_values))
      color_values <- grDevices::colorRampPalette(c("#4f9bd8", "#2e648f", "#173b5b"))(length(degree_values))
      color_values <- color_values[rank(degree_values, ties.method = "first")]
      plot(
        graph,
        layout = layout,
        vertex.size = size_values,
        vertex.color = color_values,
        vertex.frame.color = "#0b263d",
        vertex.label = label_values,
        vertex.label.color = "white",
        vertex.label.cex = 0.68,
        edge.color = grDevices::adjustcolor("#7e9bb4", alpha.f = 0.38),
        edge.width = 0.8,
        main = "STRING Protein Interaction Network — Top Connected Proteins",
        margin = 0.08
      )
      grDevices::dev.off()
      return(file.exists(output_file))
    }
  }

  plot_data <- utils::head(nodes, 20L)
  plot_data$gene_symbol <- factor(
    plot_data$gene_symbol,
    levels = rev(plot_data$gene_symbol)
  )
  graph <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = degree, y = gene_symbol)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = degree, yend = gene_symbol),
      color = "#aecbe1",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 4, color = "#2e79b7") +
    ggplot2::labs(
      title = "STRING Protein Interaction Hubs",
      subtitle = "Network view fallback: proteins ranked by interaction degree",
      x = "Interaction degree",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(output_file, graph, width = 10, height = 7, dpi = 300)
  file.exists(output_file)
}


prepare_immune_composition_plot <- function(
  result_data,
  max_cell_types = 30L,
  max_samples = 60L
) {
  result_data <- as.data.frame(result_data, stringsAsFactors = FALSE)
  if (!nrow(result_data) || ncol(result_data) < 2L) {
    return(NULL)
  }

  normalized_names <- tolower(gsub("[^a-z0-9]", "", names(result_data)))
  cell_type_candidates <- c("celltype", "cell", "immunecell", "population")
  cell_type_index <- match(cell_type_candidates, normalized_names, nomatch = 0L)
  cell_type_index <- cell_type_index[cell_type_index > 0L]
  cell_type_index <- if (length(cell_type_index)) cell_type_index[[1]] else 1L

  cell_types <- trimws(as.character(result_data[[cell_type_index]]))
  sample_indices <- setdiff(seq_len(ncol(result_data)), cell_type_index)
  sample_values <- lapply(
    result_data[sample_indices],
    function(values) suppressWarnings(as.numeric(as.character(values)))
  )
  usable_samples <- vapply(
    sample_values,
    function(values) any(is.finite(values)),
    logical(1)
  )
  sample_values <- sample_values[usable_samples]

  valid_cell_types <- !is.na(cell_types) & nzchar(cell_types)
  if (!length(sample_values) || !any(valid_cell_types)) {
    return(NULL)
  }

  abundance_matrix <- do.call(cbind, sample_values)
  abundance_matrix <- abundance_matrix[valid_cell_types, , drop = FALSE]
  rownames(abundance_matrix) <- make.unique(cell_types[valid_cell_types])
  abundance_matrix[!is.finite(abundance_matrix)] <- NA_real_

  normalized_cell_types <- tolower(gsub("[^a-z0-9]", "", rownames(abundance_matrix)))
  residual_rows <- normalized_cell_types %in% c(
    "other",
    "uncharacterizedcell",
    "uncharacterisedcell",
    "unresolved",
    "unresolvedother",
    "unknown"
  )

  residual_values <- if (any(residual_rows)) {
    colSums(abundance_matrix[residual_rows, , drop = FALSE], na.rm = TRUE)
  } else {
    numeric()
  }

  abundance_matrix <- abundance_matrix[!residual_rows, , drop = FALSE]
  if (!nrow(abundance_matrix)) return(NULL)

  cell_scores <- rowMeans(abs(abundance_matrix), na.rm = TRUE)
  cell_scores[!is.finite(cell_scores)] <- 0
  selected_cells <- names(sort(cell_scores, decreasing = TRUE))[
    seq_len(min(as.integer(max_cell_types), length(cell_scores)))
  ]
  abundance_matrix <- abundance_matrix[selected_cells, , drop = FALSE]

  total_samples <- ncol(abundance_matrix)
  sample_limit <- min(as.integer(max_samples), total_samples)
  if (total_samples > sample_limit) {
    sample_scores <- apply(abundance_matrix, 2, stats::var, na.rm = TRUE)
    sample_scores[!is.finite(sample_scores)] <- 0
    selected_positions <- sort(order(sample_scores, decreasing = TRUE)[seq_len(sample_limit)])
    selected_samples <- colnames(abundance_matrix)[selected_positions]
  } else {
    selected_samples <- colnames(abundance_matrix)
  }
  abundance_matrix <- abundance_matrix[, selected_samples, drop = FALSE]

  residual_summary <- NULL
  finite_residual <- residual_values[is.finite(residual_values)]
  if (length(finite_residual)) {
    residual_summary <- list(
      label = "Unresolved/other compartment",
      measured_columns = length(finite_residual),
      minimum = unname(min(finite_residual)),
      median = unname(stats::median(finite_residual)),
      maximum = unname(max(finite_residual))
    )
  }

  list(
    abundance_matrix = abundance_matrix,
    residual_summary = residual_summary,
    displayed_samples = ncol(abundance_matrix),
    total_samples = total_samples
  )
}


save_immune_composition_plot <- function(
  result_data,
  output_file,
  max_cell_types = 30L,
  max_samples = 60L
) {
  remove_stale_plot <- function() {
    if (file.exists(output_file)) unlink(output_file)
    FALSE
  }

  prepared <- prepare_immune_composition_plot(
    result_data = result_data,
    max_cell_types = max_cell_types,
    max_samples = max_samples
  )
  if (is.null(prepared)) return(remove_stale_plot())

  abundance_matrix <- prepared$abundance_matrix

  plot_data <- data.frame(
    Cell_Type = rep(rownames(abundance_matrix), times = ncol(abundance_matrix)),
    Sample = rep(colnames(abundance_matrix), each = nrow(abundance_matrix)),
    Abundance = as.vector(abundance_matrix),
    stringsAsFactors = FALSE
  )
  plot_data <- plot_data[is.finite(plot_data$Abundance), , drop = FALSE]
  if (!nrow(plot_data)) return(remove_stale_plot())

  plot_data$Cell_Type <- factor(
    plot_data$Cell_Type,
    levels = rev(rownames(abundance_matrix))
  )
  plot_data$Sample <- factor(
    plot_data$Sample,
    levels = colnames(abundance_matrix)
  )

  sample_note <- if (prepared$displayed_samples < prepared$total_samples) {
    paste0(
      prepared$displayed_samples,
      " most variable of ",
      prepared$total_samples,
      " samples"
    )
  } else {
    paste0("all ", prepared$total_samples, " samples")
  }

  residual_note <- if (!is.null(prepared$residual_summary)) {
    paste0(
      "Unresolved/other compartment is reported separately and excluded from the color scale: median ",
      round(100 * prepared$residual_summary$median, 1),
      "% (range ",
      round(100 * prepared$residual_summary$minimum, 1),
      "-",
      round(100 * prepared$residual_summary$maximum, 1),
      "%)."
    )
  } else {
    "No unresolved/other compartment was returned by the selected method."
  }

  graph <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Sample, y = Cell_Type, fill = Abundance)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradientn(
      colours = c("#071a2c", "#176b87", "#55c9c1", "#f4d35e"),
      labels = function(values) paste0(round(100 * values, 1), "%"),
      name = "Estimated\nfraction"
    ) +
    ggplot2::labs(
      title = "Named Immune-Cell Composition",
      subtitle = paste("quanTIseq estimated fractions across", sample_note),
      caption = paste(
        residual_note,
        "Fractions are computational estimates, not direct cell counts."
      ),
      x = "Sample",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, colour = "#4b5563"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    filename = output_file,
    plot = graph,
    width = max(9, min(16, 5 + 0.22 * ncol(abundance_matrix))),
    height = max(6, min(13, 3 + 0.30 * nrow(abundance_matrix))),
    dpi = 300,
    bg = "white"
  )

  file.exists(output_file)
}


save_drug_sensitivity_plot <- function(
  ranking,
  output_file,
  lower_is_sensitive = NA,
  max_compounds = 20L
) {
  remove_stale_plot <- function() {
    if (file.exists(output_file)) unlink(output_file)
    FALSE
  }

  ranking <- as.data.frame(ranking, stringsAsFactors = FALSE)
  required_columns <- c("Compound", "Mean_Response")
  if (!nrow(ranking) || !all(required_columns %in% names(ranking))) {
    return(remove_stale_plot())
  }

  ranking$Compound <- trimws(as.character(ranking$Compound))
  ranking$Mean_Response <- suppressWarnings(as.numeric(ranking$Mean_Response))
  valid <- nzchar(ranking$Compound) & is.finite(ranking$Mean_Response)
  ranking <- ranking[valid, , drop = FALSE]
  if (!nrow(ranking)) return(remove_stale_plot())

  if ("Rank" %in% names(ranking)) {
    ranking$Rank <- suppressWarnings(as.numeric(ranking$Rank))
    ranking <- ranking[order(ranking$Rank, na.last = TRUE), , drop = FALSE]
  } else {
    ranking <- ranking[order(ranking$Mean_Response), , drop = FALSE]
  }
  ranking <- utils::head(ranking, as.integer(max_compounds))

  plot_compound_names <- sub(
    "[[:space:]]*\\(BRD:.*\\)$",
    "",
    ranking$Compound
  )
  display_labels <- vapply(
    plot_compound_names,
    function(label) paste(strwrap(label, width = 46), collapse = "\n"),
    character(1)
  )
  ranking$Display_Compound <- factor(
    make.unique(display_labels),
    levels = rev(make.unique(display_labels))
  )

  direction_note <- if (isTRUE(lower_is_sensitive)) {
    "Lower source-metric values rank as more sensitive"
  } else if (identical(lower_is_sensitive, FALSE)) {
    "Higher source-metric values rank as more sensitive"
  } else {
    "Ranked assay response; interpret direction using the source metric"
  }

  graph <- ggplot2::ggplot(
    ranking,
    ggplot2::aes(x = Mean_Response, y = Display_Compound)
  ) +
    ggplot2::geom_col(width = 0.68, fill = "#2a9fba") +
    ggplot2::geom_vline(xintercept = 0, color = "#6b7f91", linewidth = 0.4) +
    ggplot2::labs(
      title = "Top Ranked Drug Responses",
      subtitle = direction_note,
      x = "Mean assay response",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8)
    )

  ggplot2::ggsave(
    filename = output_file,
    plot = graph,
    width = 10,
    height = 7.5,
    dpi = 300,
    bg = "white"
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

execute_wikipathways_agent <- function(
  genes,
  entrez_ids = character(),
  pvalue_cutoff = 0.05
) {
  output_directory <- file.path(
    real_agent_project_root,
    "output",
    "wikipathways"
  )
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  entrez_ids <- clean_wikipathways_entrez_ids(entrez_ids)
  if (length(entrez_ids) == 0L) {
    entrez_ids <- map_wikipathways_symbols_to_entrez(
      unique(genes)
    )
  }

  gmt_file <- fetch_wikipathways_gmt(
    file.path(output_directory, "cache")
  )
  enrichment <- calculate_wikipathways_enrichment(
    entrez_ids = entrez_ids,
    gmt_file = gmt_file,
    pvalue_cutoff = pvalue_cutoff,
    qvalue_cutoff = 0.20
  )
  results <- enrichment$results

  csv_file <- file.path(
    output_directory,
    "wikipathways_results.csv"
  )
  plot_file <- file.path(
    output_directory,
    "wikipathways_pathways.png"
  )
  readr::write_csv(results, csv_file)
  save_enrichment_dotplot(
    results,
    plot_file,
    "WikiPathways Enrichment"
  )

  interpretation_note <- if (length(unique(genes)) > 5000L) {
    paste(
      "The input contained",
      format(length(unique(genes)), big.mark = ","),
      paste(
        "genes; over-representation results are usually most meaningful",
        "for a filtered significant-gene set."
      )
    )
  } else {
    ""
  }

  list(
    success = TRUE,
    agent = "wikipathways",
    result = list(
      results = results,
      input_genes = unique(genes),
      mapped_genes = enrichment$mapped_genes
    ),
    csv = csv_file,
    plot = plot_file,
    rows = nrow(results),
    message = paste(
      "WikiPathways identified",
      nrow(results),
      "enriched pathways from",
      length(enrichment$mapped_genes),
      "mapped Entrez genes.",
      interpretation_note
    )
  )
}

execute_hallmark_agent <- function(genes, pvalue_cutoff = 0.05) {
  stopifnot(requireNamespace("msigdbr", quietly = TRUE), requireNamespace("clusterProfiler", quietly = TRUE))
  output_directory <- file.path(real_agent_project_root, "output", "hallmark")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
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
  plot_file <- file.path(output_directory, "string_network.png")
  readr::write_csv(result$top_hubs, csv_file)
  readr::write_csv(result$interactions, interaction_file)
  save_string_network_plot(result, plot_file)

  list(
    success = TRUE, agent = "string", result = result,
    csv = csv_file, interactions = interaction_file, plot = plot_file,
    rows = nrow(result$top_hubs),
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
  plot_file <- file.path(output_directory, "immune_composition_heatmap.png")
  readr::write_csv(result_table, csv_file)
  save_immune_composition_plot(result_table, plot_file)

  list(
    success = TRUE, agent = "immune", result = result,
    csv = csv_file, plot = plot_file, rows = nrow(result_table),
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
    plot_file <- file.path(output_directory, "drug_response_ranking.png")
    readr::write_csv(ranking, csv_file)
    save_drug_sensitivity_plot(ranking, plot_file, lower_is_sensitive = NA)
    return(list(success = TRUE, agent = "drug", result = list(results = ranking), csv = csv_file, plot = plot_file, rows = nrow(ranking), message = paste("Drug Sensitivity ranked", nrow(ranking), "PRISM compounds.")))
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
  plot_file <- file.path(output_directory, "drug_response_ranking.png")
  readr::write_csv(ranking, csv_file)
  save_drug_sensitivity_plot(
    ranking,
    plot_file,
    lower_is_sensitive = lower_is_sensitive
  )

  list(
    success = TRUE, agent = "drug",
    result = list(
      agent_name = "Drug Sensitivity", results = ranking,
      response_column = response_column, lower_is_sensitive = lower_is_sensitive
    ),
    csv = csv_file, plot = plot_file, rows = nrow(ranking),
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
  pvalue_cutoff = 0.05,
  effect_cutoff = 1,
  max_genes = 2000L
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

      gene_input <- prepare_gene_input(
        data,
        pvalue_cutoff = pvalue_cutoff,
        effect_cutoff = effect_cutoff,
        max_genes = max_genes
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

        wikipathways = execute_wikipathways_agent(
          genes = gene_input$genes,
          entrez_ids = unique(c(
            gene_input$entrez_ids,
            extract_wikipathways_entrez_ids(gene_input$data)
          )),
          pvalue_cutoff = pvalue_cutoff
        ),

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
      result$original_gene_count <- gene_input$original_gene_count
      result$selection_note <- gene_input$selection_note
      result$mapping_metadata <- gene_input$mapping_metadata
      result$mapping <- gene_input$mapping
      result$message <- paste(result$message, gene_input$selection_note)

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
    "hgnc_symbol",
    "gene_id",
    "geneid",
    "ensembl",
    "ensembl_gene_id",
    "entrezid"
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
    # Handle a matrix stored as samples × genes. Preserve a leading
    # sample-identifier column before numeric coercion and transpose.
    leading_values <- as.character(data[[1L]])
    leading_is_identifier <- !is.numeric(data[[1L]]) &&
      length(unique(leading_values[!is.na(leading_values) & nzchar(leading_values)])) >= 2L
    remaining_numeric <- if (ncol(data) > 1L) {
      mean(vapply(data[-1L], function(column) {
        values <- suppressWarnings(as.numeric(as.character(column)))
        mean(!is.na(values)) >= 0.8
      }, logical(1))) >= 0.8
    } else FALSE

    if (leading_is_identifier && remaining_numeric) {
      numeric_data <- data[-1L]
      numeric_data[] <- lapply(numeric_data, function(column) {
        suppressWarnings(as.numeric(as.character(column)))
      })
      expression_matrix <- t(as.matrix(numeric_data))
      colnames(expression_matrix) <- make.unique(leading_values)
      rownames(expression_matrix) <- names(numeric_data)
    } else {
      expression_matrix <- as.matrix(data)
      storage.mode(expression_matrix) <- "numeric"
      if (nrow(expression_matrix) < ncol(expression_matrix)) {
        expression_matrix <- t(expression_matrix)
      }
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

  genes <- sub("\\s*\\([0-9]+\\)\\s*$", "", genes)

  identifier_type <- detect_gene_identifier_type(genes)
  if (identifier_type %in% c("ensembl", "entrez", "mixed")) {
    mapped <- map_gene_identifiers(genes, identifier_type)
    mapping_table <- mapped$mapping_table
    keep_mapped <- !is.na(mapping_table$SYMBOL) & nzchar(trimws(mapping_table$SYMBOL))
    expression_matrix <- expression_matrix[keep_mapped, , drop = FALSE]
    genes <- mapping_table$SYMBOL[keep_mapped]
  }

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
