# Portable, self-contained reporting overrides. This module is sourced after
# results_helpers.R so legacy callers retain their function names and gain the
# production report contract without coupling reporting to Shiny UI code.

report_value <- function(value, fallback = "Not recorded") {
  if (is.null(value) || !length(value) || all(is.na(value))) return(fallback)
  value <- paste(as.character(value), collapse = ", ")
  if (!nzchar(trimws(value))) fallback else value
}

report_timestamp <- function(value) {
  if (is.null(value) || !length(value) || all(is.na(value))) return("Not recorded")
  format(as.POSIXct(value), "%Y-%m-%d %H:%M:%S %Z")
}

report_key_value_table <- function(values) {
  rows <- paste0(
    "<tr><th scope='row'>", html_escape_value(names(values)), "</th><td>",
    html_escape_value(vapply(values, report_value, character(1))), "</td></tr>",
    collapse = ""
  )
  paste0("<div class='table-wrap compact'><table><tbody>", rows, "</tbody></table></div>")
}

report_preview_table <- function(data, maximum_rows = 25L) {
  if (is.null(data) || !nrow(data)) return("<div class='empty'>No significant result rows were returned.</div>")
  preview <- format_numeric_for_display(utils::head(data, maximum_rows), digits = 2L)
  header <- paste0("<tr>", paste0("<th scope='col'>", html_escape_value(names(preview)), "</th>", collapse = ""), "</tr>")
  body <- paste0(apply(preview, 1, function(row) {
    paste0("<tr>", paste0("<td>", html_escape_value(row), "</td>", collapse = ""), "</tr>")
  }), collapse = "")
  note <- if (nrow(data) > maximum_rows) {
    paste0("<p class='note'>Previewing ", maximum_rows, " of ", format(nrow(data), big.mark = ","), " rows. The complete full-precision CSV is in the results bundle.</p>")
  } else ""
  paste0("<div class='table-wrap'><table><thead>", header, "</thead><tbody>", body, "</tbody></table></div>", note)
}

report_list_html <- function(values, ordered = FALSE, empty = "Not available for this result.") {
  values <- as.character(values %or_else% character())
  values <- values[!is.na(values) & nzchar(trimws(values))]
  if (!length(values)) return(paste0("<p class='note'>", html_escape_value(empty), "</p>"))
  tag <- if (ordered) "ol" else "ul"
  paste0("<", tag, ">", paste0("<li>", html_escape_value(values), "</li>", collapse = ""), "</", tag, ">")
}

report_deep_narrative_html <- function(value) {
  lines <- strsplit(as.character(value %or_else% ""), "\\r?\\n", perl = TRUE)[[1L]]
  lines <- lines[nzchar(trimws(lines))]
  suppress_section <- FALSE
  paste0(vapply(lines, function(line) {
    raw_text <- trimws(gsub("\\*\\*", "", line, perl = TRUE))
    known_sections <- c(
      "integrated biological interpretation", "evidence convergence and distinctions",
      "candidate regulatory and network model", "result-grounded research hypotheses",
      "testable hypotheses", "validation and next analyses", "recommended next analyses",
      "interpretation limits", "interpretive limits", "interpretive boundaries"
    )
    plain_heading <- tolower(raw_text) %in% known_sections
    heading <- grepl("^\\s*\\*\\*.*\\*\\*\\s*$", line, perl = TRUE) || plain_heading
    text <- html_escape_value(raw_text)
    if (heading) {
      suppress_section <<- tolower(text) %in% c(
        "result-grounded research hypotheses", "testable hypotheses",
        "validation and next analyses", "recommended next analyses",
        "interpretation limits", "interpretive limits", "interpretive boundaries"
      )
      if (suppress_section) return("")
    } else if (suppress_section) {
      return("")
    }
    if (heading) return(paste0("<h4>", text, "</h4>"))
    if (grepl("^[0-9]+\\.\\s+", text, perl = TRUE)) return(paste0("<div class='deep-narrative-item'>", text, "</div>"))
    if (grepl("^[*+-]\\s+", text, perl = TRUE)) return(paste0("<div class='deep-narrative-subitem'>", sub("^[*+-]\\s+", "", text, perl = TRUE), "</div>"))
    paste0("<p>", text, "</p>")
  }, character(1)), collapse = "")
}

report_plotly_library <- function() {
  path <- system.file("htmlwidgets/lib/plotlyjs/plotly-latest.min.js", package = "plotly")
  if (!nzchar(path) || !file.exists(path)) return("")
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

report_json <- function(value) {
  json <- jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", dataframe = "columns", digits = 10)
  gsub("</", "<\\/", json, fixed = TRUE)
}

report_string_3d_details_html <- function(summary, interactions = NULL) {
  ranked <- summary$data[order(summary$data$primary, decreasing = TRUE), , drop = FALSE]
  top <- utils::head(ranked, 6L)
  top_hubs <- paste0(
    "<li><span>", html_escape_value(top$label), "</span><strong>",
    html_escape_value(format(signif(top$primary, 5L), trim = TRUE)), " degree</strong></li>",
    collapse = ""
  )
  eligible_edges <- 0L
  if (!is.null(interactions) && nrow(interactions) && all(c("from_name", "to_name") %in% names(interactions))) {
    node_names <- as.character(ranked$label)
    edges <- interactions[
      as.character(interactions$from_name) %in% node_names &
        as.character(interactions$to_name) %in% node_names &
        as.character(interactions$from_name) != as.character(interactions$to_name),
      ,
      drop = FALSE
    ]
    if (nrow(edges)) {
      edge_key <- vapply(seq_len(nrow(edges)), function(index) {
        paste(sort(c(as.character(edges$from_name[[index]]), as.character(edges$to_name[[index]]))), collapse = "::")
      }, character(1))
      eligible_edges <- length(unique(edge_key))
    }
  }
  paste0(
    "<aside class='network-guide' aria-label='STRING 3D network guide'><span class='kicker'>NETWORK GUIDE</span><h4>How to read this view</h4>",
    "<div class='network-stats'><div><strong>", nrow(ranked), "</strong><span>ranked nodes considered</span></div><div><strong>", eligible_edges, "</strong><span>eligible retrieved edges</span></div></div>",
    "<ul class='network-legend'><li><i class='legend-node'></i><span><strong>Balls</strong> are submitted proteins ranked by STRING interaction degree.</span></li><li><i class='legend-edge'></i><span><strong>Lines</strong> are retrieved STRING associations among displayed proteins; they are not directional regulatory arrows.</span></li><li><i class='legend-size'></i><span><strong>Size and colour</strong> encode interaction degree within the retrieved network.</span></li></ul>",
    "<h5>Highest-degree hubs</h5><ol class='hub-rank-list'>", top_hubs, "</ol>",
    "<h5>Mouse controls</h5><p class='network-controls'><strong>Drag</strong> to rotate · <strong>scroll</strong> to zoom · <strong>hover</strong> for exact degree and edge scores · <strong>double-click</strong> to reset.</p>",
    "<p class='network-warning'><strong>Interpret carefully:</strong> connectivity prioritizes candidates; it does not prove activity, functional importance, direction, or causality.</p></aside>"
  )
}

report_interactive_chart_html <- function(key, data, title) {
  if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    return("<div class='empty'>Install plotly and jsonlite to include interactive graphs in the report.</div>")
  }
  summary <- result_visual_summary(data, key)
  if (is.null(summary)) return("<div class='empty'>No quantitative result profile was available for an interactive graph.</div>")
  chart <- summary$data
  bar_id <- paste0("plotly-bar-", key)
  colors <- if (identical(key, "drug")) {
    list(c(0, "#443983"), c(0.5, "#21918c"), c(1, "#fde725"))
  } else {
    list(c(0, "#bdd4c5"), c(0.55, "#d98468"), c(1, "#a94f3b"))
  }
  bar_trace <- list(
    type = "bar", orientation = "h", x = chart$primary, y = chart$label,
    marker = list(color = chart$secondary, colorscale = colors, showscale = TRUE, colorbar = list(title = summary$secondary_label)),
    customdata = Map(c, chart$secondary, chart$tertiary),
    hovertemplate = paste0("<b>%{y}</b><br>", summary$metric, ": %{x:.4g}<br>", summary$secondary_label, ": %{customdata[0]:.4g}<br>", summary$tertiary_label, ": %{customdata[1]:.4g}<extra></extra>")
  )
  bar_layout <- list(
    title = list(text = paste0("<b>", title, " evidence profile</b><br><sup>Hover for exact values · drag to zoom</sup>"), x = 0.02),
    xaxis = list(title = summary$metric, zeroline = FALSE, gridcolor = "#e4e9e4"),
    yaxis = list(title = "", automargin = TRUE),
    margin = list(l = 220, r = 35, t = 75, b = 60),
    paper_bgcolor = "#fffdf8", plot_bgcolor = "#fffdf8", font = list(color = "#24332d")
  )
  config <- list(displaylogo = FALSE, responsive = TRUE, scrollZoom = TRUE)
  bar_html <- paste0(
    "<div class='interactive-chart-single'><div id='", bar_id, "' class='interactive-chart'></div></div>",
    "<script>Plotly.newPlot('", bar_id, "',", report_json(list(bar_trace)), ",", report_json(bar_layout), ",", report_json(config), ");</script>"
  )
  if (!identical(key, "string")) return(bar_html)

  interactions <- if (is_existing_result_file(result_files$string$interactions)) {
    safe_result_csv(result_files$string$interactions)
  } else NULL
  three_widget <- evidence_3d_plot(summary, title, interactions)
  if (is.null(three_widget)) return(bar_html)
  three_built <- plotly::plotly_build(three_widget)
  three_id <- "plotly-3d-string"
  paste0(
    bar_html,
    "<div class='string-network-block'><h4>Connected 3D STRING interaction network</h4><div class='string-network-layout'><div id='", three_id, "' class='interactive-chart string-3d-chart'></div>",
    report_string_3d_details_html(summary, interactions), "</div></div>",
    "<script>Plotly.newPlot('", three_id, "',", report_json(three_built$x$data), ",", report_json(three_built$x$layout), ",", report_json(config), ");</script>"
  )
}

report_manifest <- function(selected_agents) {
  rows <- list()
  for (key in intersect(selected_agents, names(result_files))) {
    for (kind in intersect(c("csv", "plot", "interactions"), names(result_files[[key]]))) {
      path <- result_files[[key]][[kind]]
      exists <- is_existing_result_file(path)
      rows[[length(rows) + 1L]] <- data.frame(
        agent = key,
        artifact = kind,
        file = basename(path),
        exists = exists,
        size_bytes = if (exists) file.info(path)$size else NA_real_,
        md5 = if (exists) unname(tools::md5sum(path)) else NA_character_,
        rows = if (exists && identical(kind, "csv")) nrow(safe_result_csv(path)) else NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

report_package_versions <- function() {
  packages <- c(
    "R" = NA_character_, shiny = "shiny", clusterProfiler = "clusterProfiler",
    org.Hs.eg.db = "org.Hs.eg.db", ReactomePA = "ReactomePA", STRINGdb = "STRINGdb",
    msigdbr = "msigdbr", GSVA = "GSVA", immunedeconv = "immunedeconv",
    ggplot2 = "ggplot2", plotly = "plotly", htmlwidgets = "htmlwidgets",
    jsonlite = "jsonlite", httr2 = "httr2"
  )
  versions <- vapply(names(packages), function(label) {
    package <- packages[[label]]
    if (is.na(package)) return(paste(R.version$major, R.version$minor, sep = "."))
    if (!requireNamespace(package, quietly = TRUE)) return("Not installed")
    as.character(utils::packageVersion(package))
  }, character(1))
  data.frame(component = names(versions), version = unname(versions), stringsAsFactors = FALSE)
}

report_pathway_gene_sets <- function(keys = c("go", "kegg", "reactome", "wikipathways")) {
  sets <- list()
  for (key in intersect(keys, names(result_files))) {
    data <- safe_result_csv(result_files[[key]]$csv)
    gene_column <- find_result_column(data, c("geneID", "gene_id", "genes", "core_enrichment", "leadingEdge"))
    if (is.null(gene_column) || !nrow(data)) next
    values <- unlist(strsplit(as.character(data[[gene_column]]), "[/;,|[:space:]]+", perl = TRUE), use.names = FALSE)
    values <- unique(trimws(values[!is.na(values) & nzchar(trimws(values))]))
    if (!length(values)) next
    if (all(grepl("^[0-9]+$", values)) &&
        requireNamespace("AnnotationDbi", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      mapped <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = values,
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
      ))
      mapped <- unname(as.character(mapped))
      values <- unique(mapped[!is.na(mapped) & nzchar(mapped)])
    }
    if (length(values)) sets[[key]] <- values
  }
  sets
}

report_pathway_overlap_html <- function() {
  sets <- report_pathway_gene_sets()
  if (length(sets) < 2L) {
    return("<div class='empty'>Member-gene overlap requires at least two pathway databases with saved gene-level result columns.</div>")
  }
  all_genes <- unlist(sets, use.names = FALSE)
  frequencies <- sort(table(all_genes), decreasing = TRUE)
  shared <- frequencies[frequencies >= 2L]
  shared_text <- if (length(shared)) {
    paste0(names(utils::head(shared, 15L)), " (", as.integer(utils::head(shared, 15L)), " databases)")
  } else character()
  unique_text <- vapply(names(sets), function(key) {
    others <- unique(unlist(sets[setdiff(names(sets), key)], use.names = FALSE))
    unique_genes <- setdiff(sets[[key]], others)
    examples <- if (length(unique_genes)) paste(utils::head(unique_genes, 8L), collapse = ", ") else "none"
    paste0(agent_titles[[key]], ": ", length(unique_genes), " unique genes; examples: ", examples)
  }, character(1))
  paste0(
    "<div class='evidence-card'><p><strong>Shared member genes</strong> are genes appearing in significant result rows from two or more pathway databases.</p>",
    report_list_html(shared_text, empty = "No member gene occurred in two or more pathway result sets."),
    "<p><strong>Database-specific members</strong></p>", report_list_html(unique_text),
    "<p class='note'>This is a deterministic overlap of the complete saved significant-result tables. It is not independent replication and does not establish pathway activation.</p></div>"
  )
}

report_hub_evidence_html <- function() {
  data <- safe_result_csv(result_files$string$csv)
  symbol_column <- find_result_column(data, c("gene_symbol", "preferred_name", "symbol", "gene"))
  degree_column <- find_result_column(data, c("degree", "connectivity", "score"))
  image_uri <- image_data_uri(result_files$string$plot)
  figure <- if (!is.null(image_uri)) paste0(
    "<figure class='hub-figure'><img src='", image_uri,
    "' alt='STRING interaction network of candidate hub proteins'><figcaption>Retrieved STRING subnetwork used for connectivity ranking. Lines represent database-supported associations, not proven causal regulation.</figcaption></figure>"
  ) else "<div class='empty'>The STRING network image was not available for this run.</div>"
  rationale <- character()
  if (!is.null(symbol_column) && !is.null(degree_column) && nrow(data)) {
    degree <- suppressWarnings(as.numeric(data[[degree_column]]))
    keep <- is.finite(degree)
    ranked <- data[keep, , drop = FALSE]
    ranked$.degree <- degree[keep]
    ranked <- ranked[order(ranked$.degree, decreasing = TRUE), , drop = FALSE]
    ranked <- utils::head(ranked, 8L)
    rationale <- paste0(
      as.character(ranked[[symbol_column]]), " — interaction degree ",
      format(ranked$.degree, trim = TRUE), "; prioritized by connectivity within the retrieved network."
    )
  }
  paste0(figure, "<h4>Why these proteins are listed as hubs</h4>",
    report_list_html(rationale, empty = "Hub-rank evidence was not available."),
    "<p class='note'>Degree is a network-prioritization feature. It does not prove functional importance, direction, or causality.</p>")
}

report_chea_evidence_html <- function() {
  data <- safe_result_csv(result_files$chea$csv)
  term_column <- find_result_column(data, c("Term", "term", "regulator", "Description"))
  adjusted_column <- find_result_column(data, c("Adjusted.P.value", "p.adjust", "padj", "FDR"))
  score_column <- find_result_column(data, c("Combined.Score", "combined_score", "score"))
  if (is.null(term_column) || !nrow(data)) return("<div class='empty'>ChEA regulator evidence was not available.</div>")
  adjusted <- if (is.null(adjusted_column)) rep(NA_real_, nrow(data)) else suppressWarnings(as.numeric(data[[adjusted_column]]))
  score <- if (is.null(score_column)) rep(NA_real_, nrow(data)) else suppressWarnings(as.numeric(data[[score_column]]))
  order_index <- if (any(is.finite(adjusted))) order(adjusted, na.last = TRUE) else order(score, decreasing = TRUE, na.last = TRUE)
  selected <- utils::head(order_index, 8L)
  regulator <- sub("[[:space:]].*$", "", as.character(data[[term_column]][selected]))
  evidence <- vapply(seq_along(selected), function(index) {
    row <- selected[[index]]
    metrics <- character()
    if (is.finite(adjusted[[row]])) metrics <- c(metrics, paste0("adjusted P = ", format(adjusted[[row]], scientific = TRUE, digits = 3)))
    if (is.finite(score[[row]])) metrics <- c(metrics, paste0("combined score = ", format(score[[row]], digits = 4)))
    paste0(regulator[[index]], " — ", if (length(metrics)) paste(metrics, collapse = "; ") else "ranked ChEA target-set enrichment")
  }, character(1))
  paste0(report_list_html(evidence),
    "<p class='note'>ChEA enrichment indicates overlap with published regulator target sets. It does not demonstrate regulator activity in the submitted specimens.</p>")
}

report_consolidated_hypothesis_html <- function(synthesis) {
  convergences <- utils::head(as.character(synthesis$convergences %or_else% character()), 3L)
  hubs <- utils::head(as.character(synthesis$hub_candidates %or_else% character()), 3L)
  signal_text <- if (length(convergences)) paste(convergences, collapse = "; ") else "the recurring enriched programs"
  hub_text <- if (length(hubs)) paste(hubs, collapse = ", ") else "the highest-ranked network hubs"
  hypothesis <- paste0(
    "Working hypothesis: the recurring signals (", signal_text,
    ") define a coordinated gene program in which ", hub_text,
    " are candidate network anchors. Test this single model in an independent biological contrast using targeted perturbation or orthogonal protein/RNA measurements, with prespecified positive and negative controls."
  )
  paste0("<div class='hypothesis-card'><p>", html_escape_value(hypothesis),
    "</p><p class='note'>This is a result-grounded hypothesis for validation, not a causal or clinical conclusion.</p></div>")
}

report_interpretation_html <- function(entry) {
  if (is.null(entry)) return("<div class='empty'>Interpretation is not available.</div>")
  evidence <- entry$observed_results %or_else% entry$evidence %or_else% character()
  paste0(
    "<div class='observed-layer'><h4>Observed results · deterministic</h4><ul>",
    paste0("<li>", html_escape_value(evidence), "</li>", collapse = ""),
    "</ul></div>",
    "<div class='ai-layer'><h4>Integrated biological interpretation</h4><p>",
    html_escape_value(entry$summary %or_else% "Interpretation is not available."),
    "</p><h5>Key findings</h5>", report_list_html(entry$key_findings),
    "<h5>Biological context</h5><p>", html_escape_value(entry$biological_context %or_else% "Additional biological context is not available."),
    "<h5>Cancer relevance</h5><p>", html_escape_value(entry$cancer_relevance %or_else% "Cannot be determined from the current data."),
    "</p></div>"
  )
}

report_agent_input_label <- function(key) {
  labels <- c(
    go = "Gene list", kegg = "Gene list", reactome = "Gene list",
    wikipathways = "Gene list", string = "Gene list", hallmark = "Gene list",
    chea = "Gene list", gsva = "Expression matrix",
    immune = "Expression matrix", drug = "Drug-response table"
  )
  report_value(labels[[key]], "Scientific agent")
}

report_agent_section <- function(key, interpretation_bundle) {
  csv_path <- result_files[[key]]$csv
  data <- safe_result_csv(csv_path)
  title <- agent_titles[[key]]
  paste0(
    "<article class='agent-section' id='agent-", key, "'>",
    "<div class='agent-heading'><div><span class='kicker'>", html_escape_value(report_agent_input_label(key)), "</span><h3>", html_escape_value(title), " Analysis</h3></div><span class='row-count'>", format(nrow(data), big.mark = ","), " rows</span></div>",
    "<p class='section-copy'>", html_escape_value(agent_descriptions[[key]]), "</p>",
    "<h4>Interactive result profile</h4><p class='note'>",
    if (identical(key, "string")) "Use the horizontal bar chart for precise hub comparison and the connected 3D network to explore retrieved protein associations." else "The horizontal bar chart is the primary quantitative comparison. Hover for exact values and drag to zoom.",
    "</p>",
    report_interactive_chart_html(key, data, title),
    "<p class='note'>The complete full-precision result table is provided as a CSV in the downloadable results bundle.</p>",
    report_interpretation_html(interpretation_bundle$agents[[key]]),
    "</article>"
  )
}

report_provider_comparison_html <- function(bundle, selected_agents) {
  comparison <- bundle$comparison
  if (is.null(comparison) || is.null(comparison$ollama) || is.null(comparison$openai)) return("")
  provider_summary <- function(value, label) {
    usage <- value$usage %or_else% list()
    cost <- suppressWarnings(as.numeric(value$estimated_cost_usd %or_else% NA_real_))
    report_key_value_table(c(
      "Provider" = label,
      "State" = interpretation_display_label(value),
      "Model" = report_value(value$model),
      "Elapsed seconds" = if (is.finite(value$elapsed_seconds %or_else% NA_real_)) round(value$elapsed_seconds, 1) else "Not recorded",
      "Input tokens" = report_value(usage$input_tokens),
      "Output tokens" = report_value(usage$output_tokens),
      "Total tokens" = report_value(usage$total_tokens),
      "Estimated API cost (USD)" = if (is.finite(cost)) paste0("~$", format(cost, digits = 3, nsmall = 3)) else "Not applicable / not recorded",
      "HTTP status" = if (is.finite(value$http_status %or_else% NA_real_)) as.character(value$http_status) else "Not recorded",
      "Request ID" = report_value(value$request_id),
      "Provider diagnostic" = report_value(value$provider_error, if (value$source %in% c("ollama", "openai")) "None" else value$reason)
    ))
  }
  comparisons <- paste0(vapply(selected_agents, function(key) {
    local_entry <- comparison$ollama$agents[[key]] %or_else% list()
    remote_entry <- comparison$openai$agents[[key]] %or_else% list()
    local_summary <- if (identical(comparison$ollama$source, "ollama")) {
      local_entry$summary %or_else% "No local model narrative was available."
    } else {
      paste("Ollama did not generate model-authored prose.", comparison$ollama$reason %or_else% "")
    }
    remote_summary <- if (identical(comparison$openai$source, "openai")) {
      remote_entry$summary %or_else% "No OpenAI model narrative was available."
    } else {
      paste("OpenAI did not generate model-authored prose.", comparison$openai$reason %or_else% "")
    }
    paste0(
      "<div class='comparison-agent'><h3>", html_escape_value(agent_titles[[key]]), "</h3><div class='comparison-grid'>",
      "<div><span class='kicker'>LOCAL OLLAMA</span><p>", html_escape_value(local_summary), "</p></div>",
      "<div class='premium-comparison'><span class='kicker'>OPENAI PREMIUM</span><p>", html_escape_value(remote_summary), "</p></div>",
      "</div></div>"
    )
  }, character(1)), collapse = "")
  paste0(
    "<section id='provider-comparison'><span class='kicker'>CONTROLLED MODEL COMPARISON</span><h2>Ollama vs OpenAI Premium</h2>",
    "<p class='section-copy'>Both providers received the same structured result digest and JSON contract. Deterministic observations and named findings were rebuilt from saved result tables after generation. Prose length or fluency is not evidence of scientific correctness.</p>",
    "<div class='comparison-grid'><div><h3>Local provider</h3>", provider_summary(comparison$ollama, "Ollama"), "</div><div class='premium-comparison'><h3>Premium provider</h3>", provider_summary(comparison$openai, "OpenAI Responses API"), "</div></div>",
    "<p class='note'><strong>Primary report narrative:</strong> ", html_escape_value(toupper(comparison$primary_provider %or_else% "computed")), ". OpenAI requests transmitted only the structured result digest, used store=false, and did not upload the original input file.</p>",
    comparisons,
    "</section>"
  )
}

build_combined_html_report <- function(
  destination,
  interpretation_bundle = NULL,
  selected_agents = names(result_files),
  report_title = "OncoProfiling Multi-Agent Analysis",
  report_badge = "COMBINED REPORT",
  include_synthesis = TRUE,
  run_context = list()
) {
  selected_agents <- intersect(selected_agents, names(result_files))
  biomarker_agents <- intersect(selected_agents, setdiff(names(result_files), "drug"))
  drug_agents <- intersect(selected_agents, "drug")
  workflow <- report_value(run_context$workflow, if (length(drug_agents) && !length(biomarker_agents)) "drug" else "biomarker")
  generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  summary_data <- result_summary_rows(selected_agents)

  if (is.null(interpretation_bundle)) {
    report_data <- setNames(lapply(selected_agents, function(key) safe_result_csv(result_files[[key]]$csv)), selected_agents)
    interpretation_bundle <- generate_interpretation_bundle(report_data, settings = list(enabled = FALSE))
  }
  interpretation_bundle <- interpretation_bundle_for_report(interpretation_bundle)
  manifest <- report_manifest(selected_agents)
  packages <- report_package_versions()
  comparison_html <- report_provider_comparison_html(interpretation_bundle, selected_agents)

  summary_rows <- if (nrow(summary_data)) paste0(apply(summary_data, 1, function(row) {
    paste0("<tr><td><strong>", html_escape_value(row[["Agent"]]), "</strong></td><td>", html_escape_value(row[["Status"]]), "</td><td>", html_escape_value(row[["Rows"]]), "</td><td>", html_escape_value(row[["Plot"]]), "</td></tr>")
  }), collapse = "") else "<tr><td colspan='4'>No agents selected.</td></tr>"

  input_context <- run_context$input %or_else% list()
  import <- input_context$import %or_else% list()
  provenance <- report_key_value_table(c(
    "Workflow" = workflow,
    "Input file" = report_value(input_context$name),
    "Input MD5" = report_value(input_context$checksum_md5),
    "Input size (bytes)" = report_value(input_context$size_bytes),
    "Rows used" = report_value(input_context$rows),
    "Columns" = report_value(input_context$columns),
    "Header detected" = report_value(import$header_detected),
    "Run started" = report_timestamp(run_context$started_at),
    "Run finished" = report_timestamp(run_context$finished_at),
    "Report generated" = generated_at
  ))

  mapping <- run_context$mapping %or_else% list()
  mapping_html <- if (length(mapping)) {
    report_key_value_table(c(
      "Identifier type" = report_value(mapping$identifier_type),
      "Unique inputs" = report_value(mapping$input_count),
      "Mapped inputs" = report_value(mapping$mapped_count),
      "Unique HGNC symbols" = report_value(mapping$output_symbol_count),
      "Unmapped inputs" = report_value(mapping$unmapped_count),
      "Mapping rate" = if (is.null(mapping$mapping_rate)) "Not recorded" else paste0(round(100 * mapping$mapping_rate, 1), "%"),
      "Duplicate mappings removed" = report_value(mapping$duplicate_mappings_removed),
      "Unmapped examples" = report_value(mapping$unmapped_examples, "None")
    ))
  } else "<div class='empty'>Identifier mapping was not required or was not recorded for this workflow.</div>"

  configuration <- run_context$configuration %or_else% list()
  configuration_html <- report_key_value_table(c(
    "Experimental design / comparison" = report_value(configuration$experimental_design, "Not supplied"),
    "Detected gene column" = report_value(configuration$gene_column),
    "Detected adjusted p/FDR column" = report_value(configuration$pvalue_column),
    "Detected effect-size column" = report_value(configuration$effect_column),
    "Genes before selection" = report_value(configuration$original_gene_count),
    "Selection rule" = report_value(configuration$selection_note),
    "Gene group" = report_value(configuration$gene_direction),
    "Protein-coding only" = report_value(configuration$protein_coding_only),
    "Adjusted p/FDR cutoff" = report_value(configuration$pvalue_cutoff),
    "Absolute effect cutoff" = report_value(configuration$effect_cutoff),
    "Maximum genes" = report_value(configuration$max_genes),
    "Tested gene universe" = report_value(configuration$tested_gene_universe, "Not supplied; agent/package defaults were used.")
  ))

  synthesis <- interpretation_bundle$synthesis %or_else% list()
  synthesis_html <- if (isTRUE(include_synthesis) && length(selected_agents)) paste0(
    "<section id='integrated-interpretation'><span class='kicker'>IAN-STYLE INTEGRATED REVIEW</span><h2>",
    html_escape_value(synthesis$title %or_else% "Integrated scientific interpretation"), "</h2><p class='source'>",
    html_escape_value(interpretation_display_label(interpretation_bundle)),
    "</p><div class='executive-summary'><h3>Executive interpretation</h3><p>",
    html_escape_value(synthesis$summary %or_else% "Synthesis is not available."),
    "</p><h3>Biological and cellular interpretation</h3><p>",
    html_escape_value(synthesis$integrated_interpretation %or_else% synthesis$summary %or_else% "Not available."),
    "</p>", if (!is.null(synthesis$deep_narrative)) paste0("<h3>Detailed IAN narrative</h3><div class='deep-narrative-text'>", report_deep_narrative_html(synthesis$deep_narrative), "</div>") else "",
    "</div><div class='synthesis-grid'><div><h3>Convergent signals</h3>", report_list_html(synthesis$convergences),
    "</div><div><h3>Candidate regulatory network</h3><p>", html_escape_value(synthesis$regulatory_network %or_else% "Not available."),
    "</p><h4>ChEA evidence and prioritization rationale</h4>", report_chea_evidence_html(),
    "</div><div class='wide-card'><h3>Hub candidates and network picture</h3>", report_hub_evidence_html(),
    "</div><div><h3>Novelty and literature context</h3><p>", html_escape_value(synthesis$novelty_context %or_else% "Novelty is not verified."),
    "</p></div><div class='wide-card'><h3>Pathway member-gene overlap</h3>", report_pathway_overlap_html(),
    "</div><div class='wide-card'><h3>One consolidated result-grounded hypothesis</h3>", report_consolidated_hypothesis_html(synthesis),
    "</div></div></section>"
  ) else ""

  biomarker_html <- if (length(biomarker_agents)) paste0(
    "<section id='biomarker-results' class='workflow-section'><span class='kicker'>BIOMARKER DISCOVERY</span><h2>Biomarker evidence and interpretation</h2><p class='section-copy'>Each agent shows deterministic observed results before any interpretive context.</p>",
    paste0(vapply(biomarker_agents, report_agent_section, character(1), interpretation_bundle = interpretation_bundle), collapse = ""),
    "</section>"
  ) else ""

  drug_html <- if (length(drug_agents)) paste0(
    "<section id='drug-results' class='workflow-section drug-workflow'><span class='kicker'>DRUG SENSITIVITY</span><h2>Pharmacogenomic response</h2><p class='section-copy'>Drug-response ranking is reported independently. Assay rank does not establish mechanism, clinical benefit, or treatment suitability.</p>",
    paste0(vapply(drug_agents, report_agent_section, character(1), interpretation_bundle = interpretation_bundle), collapse = ""),
    "</section>"
  ) else ""

  package_table <- report_preview_table(packages, maximum_rows = nrow(packages))
  manifest_table <- report_preview_table(manifest, maximum_rows = max(1L, nrow(manifest)))
  interpretation_contract <- report_value(interpretation_bundle$contract_version, interpretation_contract_version)
  plotly_library <- report_plotly_library()
  responsible_html <- paste0(
    "<section id='responsible-interpretation' class='limitations'><span class='kicker'>READ BEFORE INTERPRETING RESULTS</span><h2>Responsible interpretation</h2>",
    "<ul><li>Enrichment is over-representation and does not prove pathway activation or causality.</li><li>Network centrality does not establish functional importance.</li><li>GSVA scores are relative sample-level estimates; immune deconvolution is estimated composition.</li><li>Drug-response ranks do not establish mechanism, clinical efficacy, dose, or patient suitability.</li><li>Annotation databases and installed packages are version-dependent. Unmapped genes are excluded and reported in Methods.</li><li>Literature novelty, clinical relevance, and target actionability require dedicated external evidence review.</li><li>Results require independent statistical, biological, and clinical review.</li></ul></section>"
  )

  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>", html_escape_value(report_title), "</title><style>",
    ":root{--paper:#f7f5ef;--surface:#fffdf8;--sage:#dce8de;--lavender:#e7e1f2;--ink:#24332d;--muted:#627069;--line:#d6ddd7;--accent:#c95f45;--accent-soft:#f6ded6;--drug:#6860a8}",
    "*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.58 -apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif}.page{max-width:1180px;margin:auto;padding:44px 28px 80px}",
    ".hero,section{margin:0 0 24px;padding:28px;background:var(--surface);border:1px solid var(--line);border-radius:20px;box-shadow:0 16px 42px rgba(52,67,58,.07)}.hero{background:linear-gradient(135deg,var(--sage),var(--surface) 65%)}",
    "h1{margin:10px 0 8px;font-size:38px;line-height:1.1}h2{margin:5px 0 12px;font-size:26px}h3{font-size:21px}h4{margin:24px 0 9px}h5{margin:18px 0 5px}.badge,.kicker{display:inline-block;color:var(--accent);font-size:11px;font-weight:800;letter-spacing:.12em}.badge{padding:6px 10px;background:var(--accent-soft);border-radius:999px}.muted,.note,.section-copy,.source,figcaption{color:var(--muted)}",
    ".table-wrap{overflow:auto;border:1px solid var(--line);border-radius:12px}.table-wrap.compact{max-width:900px}table{width:100%;border-collapse:collapse;background:white;font-size:13px}th,td{padding:10px 12px;border-bottom:1px solid #e8ece8;text-align:left;vertical-align:top;white-space:nowrap}thead th{background:#edf3ee}tbody th{width:240px;background:#f4f7f4}",
    ".workflow-section{border-top:7px solid #8cac93}.drug-workflow{border-top-color:var(--drug)}.agent-section{margin:24px 0;padding:24px;background:#fff;border:1px solid var(--line);border-radius:16px}.agent-heading{display:flex;align-items:flex-start;justify-content:space-between;gap:16px}.agent-heading h3{margin:4px 0}.row-count{padding:6px 9px;background:var(--sage);border-radius:999px;font-size:12px;font-weight:700}",
    "figure{margin:22px 0}img{display:block;max-width:100%;max-height:720px;margin:auto;padding:8px;background:white;border-radius:12px}figcaption{text-align:center;margin-top:8px}.empty{padding:18px;border:1px dashed #aab7ae;border-radius:12px;color:var(--muted)}",
    ".observed-layer,.ai-layer{margin-top:20px;padding:18px;border-radius:14px}.observed-layer{background:#edf5ef;border:1px solid #c9dbcd}.ai-layer{background:#f1edf7;border:1px solid #d8cfea}.observed-layer h4,.ai-layer h4{margin-top:0}.ai-layer h4{color:#51488d}",
    ".report-nav{display:flex;flex-wrap:wrap;gap:9px;margin:18px 0 0}.report-nav a{padding:8px 12px;border-radius:999px;background:#fff;color:var(--ink);text-decoration:none;border:1px solid var(--line);font-weight:700}.run-facts{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin:18px 0}.run-facts>div{padding:13px 15px;background:rgba(255,255,255,.82);border:1px solid var(--line);border-radius:12px}.run-facts span{display:block;color:var(--muted);font-size:12px}.run-facts strong{display:block;margin-top:3px;overflow-wrap:anywhere}.executive-summary{padding:20px;background:linear-gradient(135deg,#edf5ef,#f2eef7);border-radius:15px}.deep-narrative-text{line-height:1.72}.deep-narrative-text h4{margin:22px 0 7px;color:#51488d}.deep-narrative-item{margin:8px 0;padding:10px 12px;background:#fffdf8;border-left:3px solid #8cac93;border-radius:6px}.deep-narrative-subitem{margin:5px 0 5px 18px;color:var(--muted)}.synthesis-grid,.comparison-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}.synthesis-grid>div,.comparison-grid>div{padding:17px;background:#faf9f5;border:1px solid var(--line);border-radius:13px}.synthesis-grid>.wide-card{grid-column:1/-1}.premium-comparison{background:#f3edf8!important;border-color:#d4c5e4!important}.comparison-agent{margin-top:18px;padding-top:14px;border-top:1px solid var(--line)}.comparison-agent h3{margin:0 0 10px}.interactive-chart-single{max-width:980px;margin:auto}.interactive-chart{height:540px;border:1px solid var(--line);border-radius:14px;overflow:hidden}.string-network-block{margin-top:28px}.string-network-layout{display:grid;grid-template-columns:minmax(0,2fr) minmax(270px,1fr);gap:18px;align-items:stretch}.string-3d-chart{height:610px}.network-guide{padding:20px;background:linear-gradient(160deg,#edf5ef,#f7f3fb);border:1px solid #cfdcd2;border-radius:14px}.network-guide h4{margin:5px 0 15px}.network-guide h5{margin:18px 0 8px}.network-stats{display:grid;grid-template-columns:1fr 1fr;gap:8px}.network-stats>div{padding:11px;background:#fffdf8;border:1px solid var(--line);border-radius:10px}.network-stats strong{display:block;font-size:22px;color:#a94f3b}.network-stats span{display:block;color:var(--muted);font-size:11px;line-height:1.35}.network-legend{list-style:none;padding:0;margin:17px 0}.network-legend li{display:grid;grid-template-columns:22px 1fr;gap:9px;align-items:start;margin:10px 0}.network-legend i{display:block;margin-top:4px}.legend-node{width:14px;height:14px;border-radius:50%;background:#d98468;border:2px solid #fff;box-shadow:0 0 0 1px #a94f3b}.legend-edge{width:20px;height:3px;background:#6d8d7d;border-radius:4px}.legend-size{width:17px;height:17px;border-radius:50%;background:linear-gradient(135deg,#bdd4c5,#a94f3b)}.hub-rank-list{padding-left:22px}.hub-rank-list li{padding:5px 0;border-bottom:1px solid rgba(99,119,108,.16)}.hub-rank-list li span{display:inline-block;min-width:80px}.hub-rank-list li strong{font-size:12px;color:#627069}.network-controls{font-size:13px}.network-warning{padding:11px;background:#fff6ee;border-left:4px solid var(--accent);border-radius:8px;font-size:12px}.hub-figure img{max-height:680px}.evidence-card,.hypothesis-card{padding:15px;background:#fffdf8;border:1px solid var(--line);border-radius:11px}.method-callout{padding:18px;border-left:5px solid #6860a8;background:#f2eef7;border-radius:10px}.limitations{background:#fff6ee;border-left:5px solid var(--accent)}code{overflow-wrap:anywhere}@media(max-width:900px){.string-network-layout{grid-template-columns:1fr}.string-3d-chart{height:520px}}@media(max-width:800px){.synthesis-grid,.comparison-grid,.run-facts{grid-template-columns:1fr}.interactive-chart{height:430px}.string-3d-chart{height:470px}}@media(max-width:700px){.page{padding:20px 12px 48px}.hero,section{padding:20px}.agent-section{padding:16px}h1{font-size:30px}.agent-heading{flex-direction:column}th,td{font-size:11px;padding:8px}}",
    "</style>", if (nzchar(plotly_library)) paste0("<script>", plotly_library, "</script>") else "", "</head><body><main class='page'>",
    "<header class='hero'><span class='badge'>", html_escape_value(report_badge), "</span><h1>", html_escape_value(report_title), "</h1><p class='muted'>Self-contained research report · generated ", html_escape_value(generated_at), "</p>",
    "<div class='run-facts'><div><span>Input file</span><strong>", html_escape_value(report_value(input_context$name)), "</strong></div><div><span>Workflow</span><strong>", html_escape_value(workflow), "</strong></div><div><span>Mapped analysis genes</span><strong>", html_escape_value(report_value(mapping$output_symbol_count, report_value(configuration$original_gene_count))), "</strong></div></div>",
    "<h2>Run summary</h2><div class='table-wrap'><table><thead><tr><th>Agent</th><th>Status</th><th>Rows</th><th>Plot</th></tr></thead><tbody>", summary_rows, "</tbody></table></div>",
    "<nav class='report-nav' aria-label='Report sections'><a href='#integrated-interpretation'>Integrated interpretation</a>",
    if (nzchar(comparison_html)) "<a href='#provider-comparison'>Model comparison</a>" else "",
    if (length(biomarker_agents)) "<a href='#biomarker-results'>Biomarker results</a>" else "",
    if (length(drug_agents)) "<a href='#drug-results'>Drug results</a>" else "",
    "<a href='#methods'>Methods</a><a href='#responsible-interpretation'>Responsible interpretation</a></nav></header>",
    responsible_html,
    "<section><span class='kicker'>PROVENANCE</span><h2>Input and run identity</h2>", provenance, "</section>",
    "<section><span class='kicker'>INTERPRETATION STATE</span><h2>Observed results and model interpretation</h2>",
    report_key_value_table(c(
      "State" = interpretation_display_label(interpretation_bundle),
      "Provider mode" = report_value(interpretation_bundle$provider, report_value(interpretation_bundle$source)),
      "Source" = report_value(interpretation_bundle$source),
      "Model" = report_value(interpretation_bundle$model),
      "JSON contract" = interpretation_contract,
      "OpenAI request ID" = report_value(interpretation_bundle$request_id),
      "API tokens" = report_value(interpretation_bundle$usage$total_tokens),
      "Estimated API cost (USD)" = if (is.finite(interpretation_bundle$estimated_cost_usd %or_else% NA_real_)) sprintf("%.4f", interpretation_bundle$estimated_cost_usd) else "Not available"
    )),
    "<p class='note'>Observed-result bullets are computed deterministically from saved result tables. Model-generated text is a separate interpretive layer and cannot change those observations. In OpenAI modes, only the structured result digest is transmitted, the original upload is not sent, and the request uses <code>store=false</code>. Token cost is an estimate based on rates encoded in this app version; verify current platform pricing.</p></section>",
    comparison_html, synthesis_html, biomarker_html, drug_html,
    "<section id='methods'><span class='kicker'>METHODS, MAPPING, AND VERSIONS</span><h2>Reproducible methods</h2><div class='method-callout'><strong>Interpretation method.</strong> The provider-neutral prompt adapts Dr. Tyc's IAN combined-review sequence: individual-agent review, pathway integration, groundedness checking, regulator/network assessment, one consolidated hypothesis, validation planning, and high-level synthesis. Ollama and OpenAI use the same versioned JSON contract and deterministic evidence layer. Literature similarity and novelty are not asserted without a dedicated literature review.</div><h3>Gene mapping summary</h3>", mapping_html, "<h3>Input selection and experimental context</h3>", configuration_html, "<h3>Analysis environment</h3><p>Human identifier mapping uses the installed org.Hs.eg.db annotation. Enrichment statistics, multiple-testing values, and database-specific identifiers are retained in the complete CSV artifacts. Where no tested-gene universe was supplied, the relevant package or annotation collection default was used; this limits comparability across databases. Interactive horizontal bar charts are descriptive views of recorded result metrics.</p>", package_table, "</section>",
    "<section><span class='kicker'>ARTIFACT MANIFEST</span><h2>Files represented in this report</h2>", manifest_table, "</section>",
    "</main></body></html>"
  )
  writeLines(html, destination, useBytes = TRUE)
  invisible(destination)
}

create_results_bundle <- function(
  destination,
  interpretation_bundle = NULL,
  selected_agents = names(result_files),
  run_context = list()
) {
  temp_dir <- tempfile("oncoprofiling-results-")
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  selected_agents <- intersect(selected_agents, names(result_files))

  report_path <- file.path(temp_dir, "OncoProfiling_Combined_Report.html")
  build_combined_html_report(report_path, interpretation_bundle, selected_agents, run_context = run_context)
  manifest_path <- file.path(temp_dir, "artifact_manifest.csv")
  readr::write_csv(report_manifest(selected_agents), manifest_path)
  copied <- c(report_path, manifest_path)

  mapping <- run_context$mapping_table
  if (is.data.frame(mapping) && nrow(mapping)) {
    mapping_path <- file.path(temp_dir, "gene_identifier_mapping.csv")
    readr::write_csv(mapping, mapping_path)
    copied <- c(copied, mapping_path)
  }

  for (key in selected_agents) {
    for (kind in intersect(c("csv", "plot", "interactions"), names(result_files[[key]]))) {
      source_path <- result_files[[key]][[kind]]
      if (!is_existing_result_file(source_path)) next
      target_path <- file.path(temp_dir, paste0(toupper(key), "_", basename(source_path)))
      file.copy(source_path, target_path, overwrite = TRUE)
      copied <- c(copied, target_path)
    }
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(temp_dir)
  utils::zip(destination, files = basename(copied), flags = "-j")
}
