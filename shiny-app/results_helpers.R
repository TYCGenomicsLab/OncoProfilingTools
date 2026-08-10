project_root <- if (basename(getwd()) == "shiny-app") {
  normalizePath("..", mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

analysis_output_root <- file.path(project_root, "output")

dir.create(
  analysis_output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

result_files <- list(
  go = list(
    csv = file.path(analysis_output_root, "cms4", "go_results.csv"),
    plot = file.path(
      analysis_output_root,
      "visualizations",
      "go_biological_process_dotplot.png"
    )
  ),
  kegg = list(
    csv = file.path(analysis_output_root, "cms4", "kegg_results.csv"),
    plot = file.path(
      analysis_output_root,
      "visualizations",
      "kegg_pathway_dotplot.png"
    )
  ),
  gsva = list(
    csv = file.path(
      analysis_output_root,
      "gsva_bowel",
      "gsva_hallmark_scores.csv"
    ),
    plot = file.path(
      analysis_output_root,
      "gsva_bowel",
      "gsva_hallmark_heatmap.png"
    )
  ),
  chea = list(
    csv = file.path(
      analysis_output_root,
      "chea_cms4",
      "chea_results.csv"
    ),
    plot = file.path(
      analysis_output_root,
      "chea_cms4",
      "chea_tf_dotplot.png"
    )
  ),
  reactome = list(
    csv = file.path(analysis_output_root, "reactome", "reactome_results.csv"),
    plot = file.path(analysis_output_root, "reactome", "reactome_pathways.png")
  ),
  wikipathways = list(csv = file.path(analysis_output_root, "wikipathways", "wikipathways_results.csv"), plot = file.path(analysis_output_root, "wikipathways", "wikipathways_pathways.png")),
  string = list(
    csv = file.path(analysis_output_root, "string", "string_hub_proteins.csv"),
    plot = file.path(analysis_output_root, "string", "string_network.png"),
    interactions = file.path(analysis_output_root, "string", "string_interactions.csv")
  ),
  hallmark = list(csv = file.path(analysis_output_root, "hallmark", "hallmark_results.csv"), plot = file.path(analysis_output_root, "hallmark", "hallmark_pathways.png")),
  immune = list(
    csv = file.path(analysis_output_root, "immune", "immune_cell_composition.csv"),
    plot = file.path(analysis_output_root, "immune", "immune_composition_heatmap.png")
  ),
  drug = list(
    csv = file.path(analysis_output_root, "drug", "drug_sensitivity_results.csv"),
    plot = file.path(analysis_output_root, "drug", "drug_response_ranking.png")
  )
)

if (!("analysis-output" %in% names(shiny::resourcePaths()))) {
  shiny::addResourcePath(
    prefix = "analysis-output",
    directoryPath = analysis_output_root
  )
}

agent_titles <- c(
  go = "GO",
  kegg = "KEGG",
  reactome = "Reactome",
  wikipathways = "WikiPathways",
  string = "STRING",
  hallmark = "Hallmark",
  chea = "ChEA",
  gsva = "GSVA",
  immune = "Immune Deconvolution",
  drug = "Drug Sensitivity"
)

agent_descriptions <- c(
  go = "Enriched Gene Ontology biological-process terms.",
  kegg = "Significantly enriched molecular pathways.",
  reactome = "Curated pathway enrichment.",
  wikipathways = "Community pathway enrichment.",
  string = "Protein interaction network and hub ranking.",
  hallmark = "Cancer hallmark gene-set enrichment.",
  chea = "Candidate transcription-factor regulators.",
  gsva = "Sample-level pathway activity.",
  immune = "Immune-cell composition estimates.",
  drug = "Ranked compound sensitivity responses."
)

agent_groups <- list(
  gene = list(
    title = "Gene-list agents",
    description = "Seven complementary enrichment, regulation, and network views.",
    keys = c("go", "kegg", "reactome", "wikipathways", "string", "hallmark", "chea")
  ),
  expression = list(
    title = "Expression profiling",
    description = "Pathway activity and immune composition from expression matrices.",
    keys = c("gsva", "immune")
  ),
  drug = list(
    title = "Drug Sensitivity",
    description = "A dedicated pharmacogenomic response section with a pathway-exchange hook.",
    keys = "drug"
  )
)

safe_result_csv <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble())
  }

  result <- tryCatch(
    readr::read_csv(
      path,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique"
    ),
    error = function(error) {
      tibble::tibble()
    }
  )

  if (
    identical(basename(path), "chea_results.csv") &&
      nrow(result) > 0L &&
      !all(c("Term", "Adjusted.P.value", "Combined.Score") %in% names(result))
  ) {
    return(tibble::tibble())
  }

  result
}

is_existing_result_file <- function(path) {
  length(path) == 1L && !is.na(path) && nzchar(path) && file.exists(path)
}

plot_cache_token <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) return("missing")
  information <- file.info(path)
  paste(
    as.numeric(information$mtime),
    information$size,
    unname(tools::md5sum(path)),
    sep = "-"
  )
}

numeric_display_columns <- function(data) {
  names(data)[vapply(data, is.numeric, logical(1))]
}

probability_display_columns <- function(data) {
  numeric_columns <- numeric_display_columns(data)
  normalized <- tolower(gsub("[^a-z0-9]", "", numeric_columns))
  numeric_columns[grepl("^(pvalue|pvalues|padjust|adjustedpvalue|oldpvalue|oldadjustedpvalue|qvalue|fdr)$", normalized)]
}

format_numeric_for_display <- function(data, digits = 2L) {
  display <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  probability_columns <- probability_display_columns(display)
  for (column in numeric_display_columns(display)) {
    values <- display[[column]]
    display[[column]] <- ifelse(
      is.na(values),
      "",
      if (column %in% probability_columns) {
        formatC(values, format = "e", digits = digits)
      } else {
        formatC(values, format = "f", digits = digits, big.mark = ",")
      }
    )
  }
  display
}

result_datatable <- function(data) {
  table <- DT::datatable(
    data,
    rownames = FALSE,
    filter = "none",
    options = list(
      pageLength = 10,
      lengthMenu = c(5, 10, 25, 50),
      scrollX = TRUE,
      scrollY = "390px",
      scrollCollapse = TRUE,
      autoWidth = TRUE,
      dom = "<'result-table-toolbar'f>t<'result-table-footer'ip>",
      language = list(search = "Search results:", searchPlaceholder = "Pathway, term, gene…")
    ),
    class = "stripe hover compact"
  )

  numeric_columns <- numeric_display_columns(data)
  probability_columns <- probability_display_columns(data)
  decimal_columns <- setdiff(numeric_columns, probability_columns)
  if (length(decimal_columns)) {
    table <- DT::formatRound(table, columns = decimal_columns, digits = 2)
  }
  if (length(probability_columns)) {
    table <- DT::formatSignif(table, columns = probability_columns, digits = 3)
  }
  table
}


html_escape_value <- function(value) {
  htmltools::htmlEscape(as.character(if (is.null(value)) "" else value))
}

image_data_uri <- function(path) {
  if (
    length(path) != 1L || is.na(path) || !nzchar(path) ||
    !file.exists(path) || !requireNamespace("base64enc", quietly = TRUE)
  ) {
    return(NULL)
  }

  mime_type <- switch(
    tolower(tools::file_ext(path)),
    png = "image/png",
    jpg = "image/jpeg",
    jpeg = "image/jpeg",
    "application/octet-stream"
  )

  paste0(
    "data:", mime_type, ";base64,",
    base64enc::base64encode(path)
  )
}

result_summary_rows <- function(keys = names(result_files)) {
  keys <- intersect(keys, names(result_files))

  do.call(
    rbind,
    lapply(keys, function(key) {
      csv_path <- result_files[[key]]$csv
      plot_path <- result_files[[key]]$plot
      data <- safe_result_csv(csv_path)

      data.frame(
        Agent = toupper(key),
        Status = if (!file.exists(csv_path)) {
          "Not generated"
        } else if (nrow(data) == 0) {
          "Completed — 0 significant rows"
        } else {
          "Completed"
        },
        Rows = nrow(data),
        Plot = if (is_existing_result_file(plot_path)) "Available" else "Not available",
        stringsAsFactors = FALSE
      )
    })
  )
}

interpretation_entry_html <- function(entry) {
  if (is.null(entry)) return("<div class='empty'>Interpretation is not available.</div>")

  evidence <- paste0(
    "<li>", html_escape_value(entry$evidence %or_else% character()), "</li>",
    collapse = ""
  )
  limitations <- paste0(
    "<li>", html_escape_value(entry$limitations %or_else% character()), "</li>",
    collapse = ""
  )

  paste0(
    "<div class='interpretation'><h3>Biological interpretation</h3><p>",
    html_escape_value(entry$summary %or_else% "Interpretation is not available."),
    "</p><h4>Evidence from the full result set</h4><ul>", evidence,
    "</ul><h4>Cancer relevance</h4><p>",
    html_escape_value(entry$cancer_relevance %or_else% "Not assessed."),
    "</p><h4>Limitations</h4><ul>", limitations, "</ul></div>"
  )
}

build_combined_html_report <- function(
  destination,
  interpretation_bundle = NULL,
  selected_agents = names(result_files),
  report_title = "OncoProfiling Multi-Agent Analysis",
  report_badge = "COMBINED REPORT",
  include_synthesis = TRUE
) {
  selected_agents <- intersect(selected_agents, names(result_files))
  generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  summary_data <- result_summary_rows(selected_agents)

  if (is.null(interpretation_bundle)) {
    report_data <- setNames(
      lapply(selected_agents, function(key) safe_result_csv(result_files[[key]]$csv)),
      selected_agents
    )
    interpretation_bundle <- generate_interpretation_bundle(
      report_data,
      settings = list(enabled = FALSE)
    )
  }

  summary_rows <- paste0(
    apply(summary_data, 1, function(row) {
      paste0(
        "<tr><td><strong>", html_escape_value(row[["Agent"]]), "</strong></td>",
        "<td>", html_escape_value(row[["Status"]]), "</td>",
        "<td>", html_escape_value(row[["Rows"]]), "</td>",
        "<td>", html_escape_value(row[["Plot"]]), "</td></tr>"
      )
    }),
    collapse = "\n"
  )

  sections <- vapply(selected_agents, function(key) {
    csv_path <- result_files[[key]]$csv
    plot_path <- result_files[[key]]$plot
    data <- safe_result_csv(csv_path)
    title <- agent_titles[[key]]
    image_uri <- image_data_uri(plot_path)

    image_html <- if (!is.null(image_uri)) {
      paste0(
        "<figure><img src='", image_uri, "' alt='", title,
        " visualization'><figcaption>", title,
        " visualization</figcaption></figure>"
      )
    } else {
      "<div class='empty'>Visualization not available.</div>"
    }

    table_html <- if (nrow(data) > 0) {
      preview <- format_numeric_for_display(utils::head(data, 25), digits = 2L)
      header <- paste0(
        "<tr>",
        paste0("<th>", html_escape_value(names(preview)), "</th>", collapse = ""),
        "</tr>"
      )
      body <- paste0(
        apply(preview, 1, function(row) {
          paste0(
            "<tr>",
            paste0("<td>", html_escape_value(row), "</td>", collapse = ""),
            "</tr>"
          )
        }),
        collapse = "\n"
      )
      paste0(
        "<div class='table-wrap'><table><thead>", header,
        "</thead><tbody>", body, "</tbody></table></div>",
        if (nrow(data) > 25) {
          paste0("<p class='note'>Showing the first 25 of ", format(nrow(data), big.mark = ","), " rows. The complete CSV is included in the download bundle.</p>")
        } else ""
      )
    } else if (file.exists(csv_path)) {
      "<div class='empty'>The agent completed successfully but returned zero significant rows.</div>"
    } else {
      "<div class='empty'>Results have not been generated.</div>"
    }

    paste0(
      "<section><h2>", title, " Analysis</h2>",
      "<p class='section-copy'>", html_escape_value(agent_descriptions[[key]]), "</p>",
      image_html,
      "<h3>Result preview</h3>", table_html,
      interpretation_entry_html(interpretation_bundle$agents[[key]]),
      "</section>"
    )
  }, character(1))

  synthesis <- interpretation_bundle$synthesis
  synthesis_html <- paste0(
    "<section><h2>Cross-agent synthesis</h2><p class='source'>",
    html_escape_value(interpretation_display_label(interpretation_bundle)),
    "</p><p>", html_escape_value(synthesis$summary %or_else% "Synthesis is not available."), "</p>",
    "<h3>Drug–pathway exchange</h3><p>",
    html_escape_value(synthesis$drug_pathway_context %or_else% "No Drug–pathway context is available."),
    "</p></section>"
  )

  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>", html_escape_value(report_title), "</title>",
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;background:#07111f;color:#e8eef7;line-height:1.55}",
    ".page{max-width:1180px;margin:auto;padding:40px 28px 80px}",
    ".hero,section{background:#0d1b2d;border:1px solid #20324a;border-radius:18px;padding:28px;margin-bottom:24px;box-shadow:0 18px 50px rgba(0,0,0,.24)}",
    "h1{font-size:34px;margin:0 0 8px}h2{font-size:25px;margin-top:0}h3{margin-top:28px}",
    ".muted,.note,.section-copy,.source{color:#a9b8ca}.badge{display:inline-block;padding:6px 10px;border-radius:999px;background:#123f38;color:#8df1cf;font-size:13px;font-weight:700}",
    "figure{margin:22px 0}img{display:block;max-width:100%;max-height:720px;margin:auto;border-radius:12px;background:white}figcaption{text-align:center;color:#a9b8ca;margin-top:8px}",
    ".table-wrap{overflow:auto;border:1px solid #263a54;border-radius:12px}table{border-collapse:collapse;width:100%;font-size:13px;background:#0a1727}th,td{padding:10px 12px;border-bottom:1px solid #1d3048;text-align:left;vertical-align:top;white-space:nowrap}th{background:#13243a;position:sticky;top:0}.empty{padding:20px;border:1px dashed #38506e;border-radius:12px;color:#a9b8ca}",
    ".interpretation{margin-top:24px;padding:20px;border:1px solid #2c5270;border-radius:12px;background:#102840}.interpretation h3,.interpretation h4{margin-bottom:8px}.interpretation ul{padding-left:20px}",
    "</style></head><body><main class='page'>",
    "<div class='hero'><span class='badge'>", html_escape_value(report_badge), "</span>",
    "<h1>", html_escape_value(report_title), "</h1>",
    "<p class='muted'>Generated ", html_escape_value(generated_at), ". This portable report contains every selected Results Center agent, full-precision CSV files remain available separately, and displayed numeric previews use two decimals.</p>",
    "<h3>Run summary</h3><div class='table-wrap'><table><thead><tr><th>Agent</th><th>Status</th><th>Rows</th><th>Plot</th></tr></thead><tbody>", summary_rows, "</tbody></table></div></div>",
    if (isTRUE(include_synthesis)) synthesis_html else "",
    paste0(sections, collapse = "\n"),
    "</main></body></html>"
  )

  writeLines(html, destination, useBytes = TRUE)
}


build_agent_html_report <- function(destination, agent_key, interpretation_bundle = NULL) {
  stopifnot(agent_key %in% names(result_files))
  build_combined_html_report(
    destination = destination,
    interpretation_bundle = interpretation_bundle,
    selected_agents = agent_key,
    report_title = paste(agent_titles[[agent_key]], "Analysis Report"),
    report_badge = "AGENT REPORT",
    include_synthesis = FALSE
  )
}

create_results_bundle <- function(
  destination,
  interpretation_bundle = NULL,
  selected_agents = names(result_files)
) {
  temp_dir <- tempfile("oncoprofiling-results-")
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  report_path <- file.path(temp_dir, "OncoProfiling_Combined_Report.html")
  selected_agents <- intersect(selected_agents, names(result_files))
  build_combined_html_report(report_path, interpretation_bundle, selected_agents)

  copied <- c(report_path)
  for (key in selected_agents) {
    for (kind in intersect(c("csv", "plot", "interactions"), names(result_files[[key]]))) {
      source_path <- result_files[[key]][[kind]]
      if (is_existing_result_file(source_path)) {
        target_name <- paste0(toupper(key), "_", basename(source_path))
        target_path <- file.path(temp_dir, target_name)
        file.copy(source_path, target_path, overwrite = TRUE)
        copied <- c(copied, target_path)
      }
    }
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(temp_dir)
  utils::zip(destination, files = basename(copied), flags = "-j")
}

result_tab_ui <- function(key, title, description) {
  tabPanel(
    title,

    div(
      class = "result-tab-layout",

      div(
        class = "result-evidence-column",

        div(
          class = "result-visual-panel result-inner-card",

          div(
            class = "result-panel-heading",
            div(
              h3(paste(title, "Visualization")),
              p(description)
            ),
            div(
              class = "result-heading-actions",
              uiOutput(paste0(key, "_result_status")),
              actionButton(
                inputId = paste0("open_plot_", key),
                label = "Full screen",
                class = "result-fullscreen-button",
                icon = icon("expand")
              )
            )
          ),

          uiOutput(paste0(key, "_result_plot"))
        ),

        div(
          class = "result-table-panel result-inner-card",

          div(
            class = "result-panel-heading",
            div(
              h3(paste(title, "Results")),
              p("Scores use two decimals; probabilities use scientific notation. CSV files retain source precision.")
            ),

            div(
              class = "result-table-actions",
              downloadButton(
                outputId = paste0("download_report_", key),
                label = "HTML report",
                class = "result-download-button result-report-button"
              ),
              downloadButton(
                outputId = paste0("download_", key),
                label = "Download CSV",
                class = "result-download-button"
              )
            )
          ),

          DTOutput(paste0(key, "_result_table"))
        )
      ),

      div(
        class = "result-interpretation-column result-inner-card",

        div(
          class = "result-panel-heading",
          div(
            span(class = "interpretation-card-label", "LOCAL BIOLOGICAL INTERPRETATION"),
            h3(paste(title, "result-wide insight")),
            p("A concise digest computed across the complete result table—not only the first row.")
          )
        ),

        div(
          class = "biological-interpretation-card",
          uiOutput(paste0(key, "_interpretation"))
        )
      )
    )
  )
}

result_group_ui <- function(group, selected_agents) {
  keys <- intersect(group$keys, selected_agents)
  if (!length(keys)) return(NULL)

  div(
    class = paste("result-agent-group", if (identical(keys, "drug")) "result-agent-group-drug" else ""),
    div(
      class = "result-group-heading",
      div(h3(group$title), p(group$description)),
      span(paste(length(keys), if (length(keys) == 1L) "agent" else "agents"))
    ),
    do.call(
      tabsetPanel,
      c(
        list(id = paste0("results_tabs_", tolower(gsub("[^a-z]", "_", group$title))), type = "pills"),
        lapply(keys, function(key) result_tab_ui(key, agent_titles[[key]], agent_descriptions[[key]]))
      )
    )
  )
}

results_center_ui <- function(selected_agents = names(agent_titles)) {
  div(
    id = "results-center",
    class = "results-center glass-card",

    div(
      class = "results-center-heading",

      div(
        div(class = "section-label", "STEP 4"),
        h2("Analysis Results Center"),
        p(
          class = "section-description",
          paste(
            "Explore generated outputs from every selected scientific agent.",
            "Evidence stays on the left; full-result biological interpretation stays on the right."
          )
        )
      ),

      div(
        class = "results-center-actions",
        div(
          class = "results-ready-badge",
          span(class = "results-ready-dot"),
          "Saved outputs"
        ),
        downloadButton(
          "download_combined_report",
          "Combined HTML Report",
          class = "results-action-button"
        ),
        downloadButton(
          "download_all_results",
          "Download All (.zip)",
          class = "results-action-button results-action-primary"
        )
      )
    ),

    uiOutput("analysis_completion_banner"),

    div(
      class = "cross-agent-shell",
      div(
        class = "cross-agent-heading",
        div(
          span(class = "section-label", "SHARED CONTEXT"),
          h3("Cross-agent synthesis"),
          p("Selected agents publish concise, standardized outputs to one local synthesis layer.")
        )
      ),
      uiOutput("cross_agent_synthesis")
    ),

    lapply(agent_groups, result_group_ui, selected_agents = intersect(selected_agents, names(agent_titles)))
  )
}


ollama_model_display <- function(model) {
  model <- trimws(as.character(model %or_else% ""))
  if (length(model) && nzchar(model[[1L]])) model[[1L]] else "the selected model"
}

interpretation_display_label <- function(bundle) {
  model <- ollama_model_display(bundle$model)
  switch(
    bundle$source,
    loading = paste0("Loading ", model, " locally…"),
    generating = paste0("Analyzing full result table with ", model, "…"),
    pending = paste0("Analyzing full result table with ", model, "…"),
    ollama = paste("Biological interpretation generated locally with", model),
    bundle$source_label %or_else% "Rule-based fallback"
  )
}

interpretation_source_badge <- function(bundle) {
  switch(
    bundle$source,
    loading = "OLLAMA LOADING",
    generating = "OLLAMA GENERATING",
    pending = "OLLAMA GENERATING",
    ollama = "OLLAMA COMPLETE",
    "SAFE FALLBACK"
  )
}

interpretation_status_label <- function(bundle) {
  switch(
    bundle$source,
    loading = "OLLAMA LOADING",
    generating = "OLLAMA GENERATING",
    pending = "OLLAMA GENERATING",
    ollama = "OLLAMA COMPLETE",
    "Rule summary active"
  )
}

build_ollama_progress_bundle <- function(exchanges, state, model) {
  stopifnot(state %in% c("loading", "generating"))
  reason <- switch(
    state,
    loading = paste(
      "Ollama is starting locally.",
      "Scientific results and downloads remain available."
    ),
    generating = paste(
      "Ollama is generating a private, result-wide interpretation locally.",
      "Scientific results and downloads remain available."
    )
  )
  bundle <- build_rule_interpretation_bundle(exchanges, reason)
  bundle$source <- state
  bundle$model <- ollama_model_display(model)
  bundle$source_label <- interpretation_display_label(bundle)
  bundle
}

interpretation_source_ui <- function(bundle) {
  local_model <- bundle$source %in% c("loading", "generating", "pending", "ollama")
  div(
    class = paste("interpretation-source", paste0("interpretation-source-", bundle$source)),
    span(interpretation_source_badge(bundle)),
    strong(interpretation_display_label(bundle)),
    if (!is.null(bundle$model)) {
      tags$small(if (local_model) paste(bundle$model, "· Running locally") else bundle$model)
    },
    p(bundle$reason)
  )
}

build_agent_interpretation <- function(agent_id, data, entry = NULL, bundle = NULL) {
  if (is.null(bundle) || is.null(entry)) {
    fallback <- generate_interpretation_bundle(
      setNames(list(data), agent_id),
      settings = list(enabled = FALSE)
    )
    bundle <- fallback
    entry <- fallback$agents[[agent_id]]
  }

  if (is.null(data) || nrow(data) == 0L) {
    return(
      tagList(
        interpretation_source_ui(bundle),
        tags$p(
          class = "interpretation-empty",
          paste(
            "No significant result rows are available for this agent.",
            "This can be a valid biological outcome depending on the input data and thresholds."
          )
        )
      )
    )
  }

  evidence_items <- lapply(entry$evidence %or_else% character(), tags$li)
  limitation_items <- lapply(entry$limitations %or_else% character(), tags$li)

  tagList(
    interpretation_source_ui(bundle),
    div(
      class = "interpretation-section",
      h4("Biological summary"),
      p(entry$summary)
    ),
    tags$details(
      class = "interpretation-details",
      tags$summary("Full-result evidence"),
      div(
      class = "interpretation-section",
      h4("Coverage and distributions"),
      tags$ul(evidence_items)
      )
    ),
    div(
      class = "interpretation-section",
      h4("Cancer relevance"),
      p(entry$cancer_relevance)
    ),
    tags$details(
      class = "interpretation-details interpretation-limitations",
      tags$summary("Confidence and limitations"),
      div(class = "interpretation-section", tags$ul(limitation_items))
    )
  )
}

build_cross_agent_synthesis_ui <- function(bundle) {
  synthesis <- bundle$synthesis
  convergences <- lapply(synthesis$convergences %or_else% character(), tags$li)
  limitations <- lapply(synthesis$limitations %or_else% character(), tags$li)
  bridge_available <- isTRUE(synthesis$bridge$available)

  div(
    class = "cross-agent-content",
    interpretation_source_ui(bundle),
    div(class = "cross-agent-summary", p(synthesis$summary)),
    div(
      class = "cross-agent-grid",
      div(h4("Convergent signals"), tags$ul(convergences)),
      div(
        class = if (bridge_available) "drug-pathway-bridge bridge-ready" else "drug-pathway-bridge",
        h4("Drug ↔ pathway exchange"),
        p(synthesis$drug_pathway_context)
      ),
      div(h4("Interpretation limits"), tags$ul(limitations))
    )
  )
}

start_local_interpretation_job <- function(data_by_agent, settings) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("Background interpretation requires the processx package.", call. = FALSE)
  }

  worker_path <- normalizePath("run_interpretation_worker.R", mustWork = TRUE)
  helper_path <- normalizePath("interpretation.R", mustWork = TRUE)
  input_path <- tempfile("oncoprofiling-interpretation-input-", fileext = ".rds")
  output_path <- tempfile("oncoprofiling-interpretation-output-", fileext = ".rds")
  stdout_path <- tempfile("oncoprofiling-interpretation-stdout-", fileext = ".log")
  stderr_path <- tempfile("oncoprofiling-interpretation-stderr-", fileext = ".log")

  saveRDS(list(data_by_agent = data_by_agent, settings = settings), input_path)
  process <- processx::process$new(
    command = file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", worker_path, input_path, output_path, helper_path),
    stdout = stdout_path,
    stderr = stderr_path,
    cleanup = TRUE,
    cleanup_tree = TRUE
  )

  list(
    process = process,
    input_path = input_path,
    output_path = output_path,
    stdout_path = stdout_path,
    stderr_path = stderr_path
  )
}

cleanup_local_interpretation_job <- function(job, terminate = FALSE) {
  if (is.null(job)) return(invisible(NULL))
  if (isTRUE(terminate) && !is.null(job$process)) {
    try(if (isTRUE(job$process$is_alive())) job$process$kill_tree(), silent = TRUE)
  }
  paths <- unlist(job[c("input_path", "output_path", "stdout_path", "stderr_path")], use.names = FALSE)
  paths <- paths[!is.na(paths) & nzchar(paths) & file.exists(paths)]
  if (length(paths)) unlink(paths, force = TRUE)
  invisible(NULL)
}


register_results_server <- function(
  input,
  output,
  session,
  active_agents = NULL,
  analysis_complete = NULL,
  ollama_settings = NULL
) {

  active_keys <- shiny::reactive({
    keys <- if (is.null(active_agents)) names(agent_titles) else active_agents()
    intersect(keys, names(agent_titles))
  })

  result_watch_paths <- unique(unlist(
    lapply(
      result_files,
      function(agent_files) {
        unname(unlist(agent_files, use.names = FALSE))
      }
    ),
    use.names = FALSE
  ))
  result_watch_paths <- result_watch_paths[!is.na(result_watch_paths) & nzchar(result_watch_paths)]

  result_signature <- shiny::reactivePoll(
    intervalMillis = 750,
    session = session,

    checkFunc = function() {
      information <- file.info(result_watch_paths)

      paste(
        file.exists(result_watch_paths),
        information$size,
        as.numeric(information$mtime),
        collapse = "|"
      )
    },

    valueFunc = function() {
      Sys.time()
    }
  )

  interpretation_cache <- new.env(parent = emptyenv())
  interpretation_cache$key <- NULL
  interpretation_cache$value <- NULL
  interpretation_job <- new.env(parent = emptyenv())
  interpretation_job$value <- NULL
  interpretation_version <- shiny::reactiveVal(0L)

  publish_interpretation <- function(bundle, cache_key) {
    interpretation_cache$key <- cache_key
    interpretation_cache$value <- bundle
    interpretation_version(shiny::isolate(interpretation_version()) + 1L)
    invisible(bundle)
  }

  shiny::observe({
    result_signature()
    keys <- active_keys()

    if (!length(keys)) {
      cleanup_local_interpretation_job(interpretation_job$value, terminate = TRUE)
      interpretation_job$value <- NULL
      publish_interpretation(
        build_rule_interpretation_bundle(list(), "No agents are selected for synthesis."),
        "no-agents"
      )
      return()
    }

    settings <- if (is.null(ollama_settings)) {
      default_ollama_settings()
    } else {
      ollama_settings()
    }
    settings <- normalise_ollama_settings(settings)

    run_is_complete <- if (is.null(analysis_complete)) TRUE else isTRUE(analysis_complete())
    csv_paths <- vapply(keys, function(key) result_files[[key]]$csv, character(1))
    information <- file.info(csv_paths)
    cache_key <- paste(
      keys,
      file.exists(csv_paths),
      information$size,
      as.numeric(information$mtime),
      settings$enabled,
      settings$host,
      settings$model,
      settings$timeout_seconds,
      settings$num_predict,
      run_is_complete,
      collapse = "|"
    )

    if (identical(cache_key, interpretation_cache$key)) {
      return()
    }

    cleanup_local_interpretation_job(interpretation_job$value, terminate = TRUE)
    interpretation_job$value <- NULL

    data_by_agent <- setNames(
      lapply(keys, function(key) safe_result_csv(result_files[[key]]$csv)),
      keys
    )

    exchanges <- Map(build_agent_exchange, names(data_by_agent), data_by_agent)
    if (!run_is_complete) {
      publish_interpretation(
        build_rule_interpretation_bundle(
          exchanges,
          "Scientific workers are still running. Local interpretation starts automatically when the selected analyses finish."
        ),
        cache_key
      )
      return()
    }

    if (!isTRUE(settings$enabled)) {
      publish_interpretation(
        build_rule_interpretation_bundle(
          exchanges,
          "Local Ollama is disabled; the deterministic full-result summary is shown."
        ),
        cache_key
      )
      return()
    }

    validation_error <- validate_ollama_settings(settings)
    if (!is.null(validation_error)) {
      publish_interpretation(
        build_rule_interpretation_bundle(exchanges, paste("Local Ollama was not started:", validation_error)),
        cache_key
      )
      return()
    }

    loading <- build_ollama_progress_bundle(
      exchanges,
      state = "loading",
      model = settings$model
    )
    publish_interpretation(loading, cache_key)

    interpretation_job$value <- tryCatch(
      {
        job <- start_local_interpretation_job(data_by_agent, settings)
        job$key <- cache_key
        job
      },
      error = function(error) {
        fallback <- build_rule_interpretation_bundle(exchanges, friendly_ollama_error(error))
        publish_interpretation(fallback, cache_key)
        NULL
      }
    )
  })

  finalize_interpretation_job <- function() {
    job <- interpretation_job$value
    if (is.null(job)) return(invisible(FALSE))

    is_alive <- isTRUE(tryCatch(job$process$is_alive(), error = function(error) FALSE))
    if (is_alive) {
      current <- interpretation_cache$value
      if (
        identical(job$key, interpretation_cache$key) &&
        identical(current$source, "loading")
      ) {
        generating <- build_ollama_progress_bundle(
          current$exchanges %or_else% list(),
          state = "generating",
          model = current$model
        )
        publish_interpretation(generating, job$key)
      }
      return(invisible(FALSE))
    }

    bundle <- if (file.exists(job$output_path)) {
      tryCatch(readRDS(job$output_path), error = function(error) NULL)
    } else {
      NULL
    }

    if (is.null(bundle)) {
      current <- interpretation_cache$value
      exchanges <- current$exchanges %or_else% list()
      bundle <- build_rule_interpretation_bundle(
        exchanges,
        "The background local-model process ended without a valid response. Scientific results remain available with the safe rule-based summary."
      )
    }

    completed_key <- job$key
    cleanup_local_interpretation_job(job)
    interpretation_job$value <- NULL
    if (identical(completed_key, interpretation_cache$key)) {
      publish_interpretation(bundle, completed_key)
    }
    invisible(TRUE)
  }

  shiny::observe({
    shiny::invalidateLater(500, session)
    finalize_interpretation_job()
  })

  session$onSessionEnded(function() {
    cleanup_local_interpretation_job(interpretation_job$value, terminate = TRUE)
  })

  interpretation_bundle <- shiny::reactive({
    interpretation_version()
    interpretation_cache$value %or_else%
      build_rule_interpretation_bundle(list(), "No result summaries are available yet.")
  })


  output$analysis_completion_banner <- renderUI({
    result_signature()
    keys <- active_keys()
    all_summary <- result_summary_rows()
    summary_data <- all_summary[match(toupper(keys), all_summary$Agent), , drop = FALSE]
    summary_data <- summary_data[!is.na(summary_data$Agent), , drop = FALSE]
    completed <- sum(summary_data$Status != "Not generated")
    total_rows <- sum(summary_data$Rows)
    nonempty <- sum(summary_data$Rows > 0L)
    empty <- max(0L, completed - nonempty)
    bundle <- interpretation_bundle()
    interpretation_status <- interpretation_status_label(bundle)

    div(
      class = paste(
        "analysis-completion-banner",
        if (completed == length(keys)) "analysis-completion-complete" else "analysis-completion-partial"
      ),
      div(
        class = "analysis-completion-icon",
        if (completed == length(keys)) "✓" else "↻"
      ),
      div(
        class = "completion-status-copy",
        h4(if (completed == length(keys)) "Selected analysis outputs are ready" else "Analysis outputs are partially available"),
        p(
          paste0(
            completed, " of ", length(keys), " selected agents generated output files · ",
            format(total_rows, big.mark = ","), " total result rows",
            if (empty > 0L) paste0(" · ", empty, " completed with no significant terms") else ""
          )
        )
      ),
      span(
        class = paste("interpretation-status-chip", paste0("interpretation-status-", bundle$source)),
        interpretation_status
      )
    )
  })

  output$cross_agent_synthesis <- renderUI({
    build_cross_agent_synthesis_ui(interpretation_bundle())
  })

  output$download_combined_report <- downloadHandler(
    filename = function() {
      paste0("OncoProfiling_Combined_Report_", format(Sys.Date(), "%Y%m%d"), ".html")
    },
    content = function(file) {
      finalize_interpretation_job()
      current <- interpretation_cache$value %or_else% interpretation_bundle()
      build_combined_html_report(file, current, active_keys())
    },
    contentType = "text/html"
  )

  output$download_all_results <- downloadHandler(
    filename = function() {
      paste0("OncoProfiling_All_Results_", format(Sys.Date(), "%Y%m%d"), ".zip")
    },
    content = function(file) {
      finalize_interpretation_job()
      current <- interpretation_cache$value %or_else% interpretation_bundle()
      create_results_bundle(file, current, active_keys())
    },
    contentType = "application/zip"
  )

  for (agent_key in names(agent_titles)) {
    local({
      key <- agent_key
      title <- agent_titles[[key]]
      csv_path <- result_files[[key]]$csv
      plot_path <- result_files[[key]]$plot

      observeEvent(input[[paste0("open_plot_", key)]], {
        if (!is_existing_result_file(plot_path)) {
          showNotification(paste(title, "visualization is not available yet."), type = "warning")
          return()
        }

        relative_plot_path <- sub(
          paste0("^", analysis_output_root, "/?"),
          "",
          plot_path
        )

        plot_source <- image_data_uri(plot_path)
        if (is.null(plot_source)) {
          plot_source <- paste0(
            "analysis-output/", relative_plot_path,
            "?v=", plot_cache_token(plot_path)
          )
        }

        showModal(
          modalDialog(
            title = paste(title, "Visualization"),
            tags$img(
              src = plot_source,
              alt = paste(title, "analysis visualization"),
              class = "result-modal-image"
            ),
            size = "l",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      }, ignoreInit = TRUE)

      output[[paste0(key, "_result_status")]] <- renderUI({
    result_signature()
        result_data <- safe_result_csv(csv_path)

        if (!file.exists(csv_path)) {
          return(
            span(
              class = "result-status result-status-missing",
              "File missing"
            )
          )
        }

        if (nrow(result_data) == 0) {
          return(
            span(
              class = "result-status result-status-empty",
              "0 significant results"
            )
          )
        }

        span(
          class = "result-status result-status-ready",
          paste(
            format(nrow(result_data), big.mark = ","),
            "rows"
          )
        )
      })

      output[[paste0(key, "_interpretation")]] <- renderUI({
        result_signature()
        result_data <- safe_result_csv(csv_path)
        bundle <- interpretation_bundle()

        build_agent_interpretation(
          agent_id = key,
          data = result_data,
          entry = bundle$agents[[key]],
          bundle = bundle
        )
      })

      output[[paste0(key, "_result_plot")]] <- renderUI({
    result_signature()
        if (!is_existing_result_file(plot_path)) {
          result_data <- safe_result_csv(csv_path)
          message_text <- if (file.exists(csv_path) && nrow(result_data) == 0L) {
            paste(
              title, "completed successfully, but no statistically significant",
              "results passed the current threshold, so a plot would be misleading."
            )
          } else if (key == "string" && file.exists(csv_path)) {
            "STRING returned hub data, but a connected subnetwork could not be rendered. The table and edge-list download remain available."
          } else {
            paste(title, "visualization is not available yet.")
          }

          return(
            div(
              class = "result-empty-state",
              div(class = "result-empty-icon", "⌁"),
              h4("No visualization available"),
              p(message_text)
            )
          )
        }

        relative_plot_path <- sub(
          paste0("^", analysis_output_root, "/?"),
          "",
          plot_path
        )

        plot_source <- image_data_uri(plot_path)
        if (is.null(plot_source)) {
          plot_source <- paste0(
            "analysis-output/", relative_plot_path,
            "?v=", plot_cache_token(plot_path)
          )
        }

        tags$img(
          src = plot_source,
          alt = paste(title, "analysis visualization"),
          class = "result-plot-image"
        )
      })

      output[[paste0(key, "_result_table")]] <- renderDT({
    result_signature()
        result_data <- safe_result_csv(csv_path)

        validate(
          need(
            file.exists(csv_path),
            paste(title, "result file was not found.")
          ),
          need(
            nrow(result_data) > 0,
            paste(
              title,
              "completed but returned no significant result rows."
            )
          )
        )

        result_datatable(result_data)
      })

      output[[paste0("download_", key)]] <- downloadHandler(
        filename = function() {
          paste0(
            "OncoProfiling_",
            toupper(key),
            "_results.csv"
          )
        },

        content = function(file) {
          validate(
            need(
              file.exists(csv_path),
              paste(title, "result file is unavailable.")
            )
          )

          file.copy(
            from = csv_path,
            to = file,
            overwrite = TRUE
          )
        }
      )

      output[[paste0("download_report_", key)]] <- downloadHandler(
        filename = function() {
          paste0("OncoProfiling_", toupper(key), "_Report_", format(Sys.Date(), "%Y%m%d"), ".html")
        },
        content = function(file) {
          finalize_interpretation_job()
          current <- interpretation_cache$value %or_else% interpretation_bundle()
          build_agent_html_report(file, key, current)
        },
        contentType = "text/html"
      )
    })
  }

  list(interpretation_bundle = interpretation_bundle)
}
