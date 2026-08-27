project_root <- if (basename(getwd()) == "shiny-app") {
  normalizePath("..", mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

analysis_output_root <- file.path(project_root, "output")
completed_interpretation_cache_path <- file.path(
  analysis_output_root,
  "interpretation",
  "completed_interpretation.rds"
)

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
    "</p><h4>What the pattern may mean biologically</h4><p>",
    html_escape_value(entry$biological_context %or_else% "Additional biological context is not available."),
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
  interpretation_bundle <- interpretation_bundle_for_report(interpretation_bundle)

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
  visualization_copy <- if (identical(key, "string")) {
    paste(description, "Compare hub degree in the bar view, then inspect retrieved protein associations in the connected 3D network.")
  } else {
    paste(description, "Use the professional bar view for clear, exact result comparison.")
  }

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
              p(visualization_copy)
            ),
            div(
              class = "result-heading-actions",
              uiOutput(paste0(key, "_result_status")),
              downloadButton(
                outputId = paste0("download_report_", key),
                label = "HTML report",
                class = "result-download-button result-report-button"
              ),
              downloadButton(
                outputId = paste0("download_", key),
                label = "Download CSV",
                class = "result-download-button"
              ),
              actionButton(
                inputId = paste0("open_plot_", key),
                label = "Full screen",
                class = "result-fullscreen-button",
                icon = icon("expand")
              )
            )
          ),
          uiOutput(paste0(key, "_result_plot"))
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

results_center_ui <- function(selected_agents = names(agent_titles), workflow = "biomarker") {
  biomarker_workflow <- !identical(workflow, "drug")
  div(
    id = "results-center",
    class = "results-center glass-card",

    div(
      class = "results-center-heading",

      div(
        div(class = "section-label", if (biomarker_workflow) "BIOMARKER RESULTS" else "DRUG-SENSITIVITY RESULTS"),
        h2(if (biomarker_workflow) "Biomarker Results Center" else "Drug Sensitivity Results Center"),
        p(
          class = "section-description",
          paste(
            "Explore generated outputs from every selected scientific agent.",
            "Observed results stay distinct from the local interpretive layer."
          )
        )
      ),

      div(
        class = "results-center-actions",
        actionLink(
          "results_back_home",
          "← Main menu",
          class = "results-action-button results-navigation-button"
        ),
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

    if (length(selected_agents) >= 1L) div(
      class = "cross-agent-shell",
      div(
        class = "cross-agent-heading",
        div(
          span(class = "section-label", if (length(selected_agents) > 1L) "SHARED CONTEXT" else "INTERPRETATION OVERVIEW"),
          h3(if (length(selected_agents) > 1L) "Cross-agent synthesis" else "Deep interpretation overview"),
          p(if (length(selected_agents) > 1L) "Selected agents publish standardized outputs to one grounded local synthesis layer." else "The local model develops the agent result into a structured scientific narrative and validation plan.")
        )
      ),
      uiOutput("cross_agent_synthesis")
    ) else NULL,

    lapply(agent_groups, result_group_ui, selected_agents = intersect(selected_agents, names(agent_titles)))
  )
}


ollama_model_display <- function(model) {
  model <- trimws(as.character(model %or_else% ""))
  if (length(model) && nzchar(model[[1L]])) model[[1L]] else "the selected model"
}

interpretation_display_label <- function(bundle) {
  model <- ollama_model_display(bundle$model)
  provider <- bundle$provider %or_else% "ollama"
  if (identical(bundle$source, "loading")) {
    return(if (identical(provider, "ollama")) paste0("Loading ", model, " locally…") else paste0("Preparing ", model, "…"))
  }
  if (bundle$source %in% c("generating", "pending")) {
    return(if (identical(provider, "ollama")) paste0("Analyzing full result table with ", model, "…") else paste0("Analyzing the structured result digest with ", model, "…"))
  }
  switch(
    bundle$source,
    ollama = paste("Biological interpretation generated locally with", model),
    openai = paste("Premium biological interpretation generated with", model),
    comparison = paste("Ollama and OpenAI comparison completed ·", model),
    ollama_timeout = "Local interpretation timed out · computed summary available",
    ollama_unavailable = "Local Ollama unavailable · computed summary available",
    ollama_error = "Local interpretation ended · computed summary available",
    openai_timeout = "OpenAI Premium timed out · computed summary available",
    openai_unavailable = "OpenAI Premium unavailable · computed summary available",
    openai_error = "OpenAI Premium ended · computed summary available",
    comparison_timeout = "Provider comparison timed out · computed summary available",
    comparison_error = "Provider comparison ended · computed summary available",
    bundle$source_label %or_else% "Computed scientific summary"
  )
}

interpretation_source_badge <- function(bundle) {
  if (identical(bundle$source, "loading")) {
    return(switch(bundle$provider %or_else% "ollama", openai = "OPENAI PREMIUM LOADING", compare = "COMPARISON LOADING", "OLLAMA LOADING"))
  }
  if (bundle$source %in% c("generating", "pending")) {
    return(switch(bundle$provider %or_else% "ollama", openai = "OPENAI PREMIUM GENERATING", compare = "COMPARISON GENERATING", "OLLAMA GENERATING"))
  }
  switch(
    bundle$source,
    loading = "OLLAMA LOADING",
    generating = "OLLAMA GENERATING",
    pending = "OLLAMA GENERATING",
    ollama = "OLLAMA COMPLETE",
    openai = "OPENAI PREMIUM COMPLETE",
    comparison = "COMPARISON COMPLETE",
    ollama_timeout = "OLLAMA TIMED OUT",
    ollama_unavailable = "OLLAMA UNAVAILABLE",
    ollama_error = "OLLAMA ERROR",
    openai_timeout = "OPENAI TIMED OUT",
    openai_unavailable = "OPENAI UNAVAILABLE",
    openai_error = "OPENAI ERROR",
    comparison_timeout = "COMPARISON TIMED OUT",
    comparison_error = "COMPARISON ERROR",
    "COMPUTED SUMMARY"
  )
}

interpretation_status_label <- function(bundle) {
  if (identical(bundle$source, "loading")) {
    return(switch(bundle$provider %or_else% "ollama", openai = "OPENAI PREMIUM LOADING", compare = "COMPARISON LOADING", "OLLAMA LOADING"))
  }
  if (bundle$source %in% c("generating", "pending")) {
    return(switch(bundle$provider %or_else% "ollama", openai = "OPENAI PREMIUM GENERATING", compare = "COMPARISON GENERATING", "OLLAMA GENERATING"))
  }
  switch(
    bundle$source,
    loading = "OLLAMA LOADING",
    generating = "OLLAMA GENERATING",
    pending = "OLLAMA GENERATING",
    ollama = "OLLAMA COMPLETE",
    openai = "OPENAI PREMIUM COMPLETE",
    comparison = "COMPARISON COMPLETE",
    ollama_timeout = "OLLAMA TIMED OUT",
    ollama_unavailable = "OLLAMA UNAVAILABLE",
    ollama_error = "OLLAMA ERROR",
    openai_timeout = "OPENAI TIMED OUT",
    openai_unavailable = "OPENAI UNAVAILABLE",
    openai_error = "OPENAI ERROR",
    comparison_timeout = "COMPARISON TIMED OUT",
    comparison_error = "COMPARISON ERROR",
    "COMPUTED SUMMARY"
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

build_provider_progress_bundle <- function(exchanges, state, settings) {
  settings <- normalise_interpretation_settings(settings)
  models <- switch(
    settings$provider,
    openai = settings$openai_model,
    compare = paste(settings$model, "vs", settings$openai_model),
    settings$model
  )
  reason <- switch(
    state,
    loading = if (identical(settings$provider, "ollama")) {
      "Ollama is starting locally. Scientific results and downloads remain available."
    } else {
      "The selected interpretation provider is starting. Only the structured scientific result digest will be processed; scientific results remain available."
    },
    generating = if (identical(settings$provider, "ollama")) {
      "Ollama is generating a private result-wide interpretation locally. Scientific results remain available."
    } else if (identical(settings$provider, "openai")) {
      "OpenAI Premium is processing the structured result digest externally with store=false. Raw uploaded files are not transmitted."
    } else {
      "Ollama and OpenAI Premium are processing the same structured digest for a controlled comparison. Raw uploaded files are not transmitted."
    }
  )
  bundle <- build_rule_interpretation_bundle(exchanges, reason)
  bundle$source <- state
  bundle$model <- models
  bundle$provider <- settings$provider
  bundle$source_label <- interpretation_display_label(bundle)
  bundle
}

is_interpretation_progress_bundle <- function(bundle) {
  is.list(bundle) && bundle$source %in% c("loading", "generating", "pending")
}

is_terminal_interpretation_bundle <- function(bundle) {
  is.list(bundle) &&
    !is_interpretation_progress_bundle(bundle) &&
    is.list(bundle$agents) &&
    is.list(bundle$synthesis)
}

persist_completed_interpretation <- function(
  cache_key,
  bundle,
  path = completed_interpretation_cache_path
) {
  if (!is_terminal_interpretation_bundle(bundle)) return(invisible(FALSE))

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    "completed-interpretation-",
    tmpdir = dirname(path),
    fileext = ".rds.tmp"
  )
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path, force = TRUE), add = TRUE)

  saveRDS(
    list(version = 3L, cache_key = cache_key, bundle = bundle),
    temporary_path
  )
  published <- file.rename(temporary_path, path)
  if (!published) {
    published <- file.copy(temporary_path, path, overwrite = TRUE)
  }
  invisible(isTRUE(published))
}

read_completed_interpretation <- function(
  cache_key,
  path = completed_interpretation_cache_path
) {
  if (!file.exists(path)) return(NULL)
  cached <- tryCatch(readRDS(path), error = function(error) NULL)
  if (
    !is.list(cached) ||
      !identical(cached$version, 3L) ||
      !identical(cached$cache_key, cache_key) ||
      !is_terminal_interpretation_bundle(cached$bundle)
  ) {
    return(NULL)
  }
  cached$bundle
}

interpretation_bundle_for_report <- function(bundle) {
  if (!is_interpretation_progress_bundle(bundle)) return(bundle)
  build_rule_interpretation_bundle(
    bundle$exchanges %or_else% list(),
    paste(
      "The report uses the complete deterministic interpretation because",
      "model generation had not reached a terminal state when it was exported. No provisional AI text is included."
    )
  )
}

interpretation_source_ui <- function(bundle) {
  running <- bundle$source %in% c("loading", "generating", "pending")
  local_model <- bundle$source %in% c(
    "ollama",
    "ollama_timeout",
    "ollama_unavailable",
    "ollama_error"
  )
  shiny::div(
    class = paste("interpretation-source", paste0("interpretation-source-", bundle$source)),
    shiny::span(interpretation_source_badge(bundle)),
    shiny::strong(interpretation_display_label(bundle)),
    if (!is.null(bundle$model)) {
      shiny::tags$small(if (running) {
        paste(bundle$model, if (identical(bundle$provider, "ollama")) "· Running locally" else "· Processing selected provider mode")
      } else if (local_model) {
        paste(
          bundle$model,
          if (identical(bundle$source, "ollama")) "· Generated locally" else "· Local attempt ended"
        )
      } else if (bundle$source %in% c("openai", "comparison")) {
        paste(bundle$model, if (identical(bundle$source, "openai")) "· External API · store=false" else "· Controlled provider comparison")
      } else bundle$model)
    },
    shiny::p(bundle$reason)
  )
}

provider_comparison_agent_ui <- function(bundle, agent_id) {
  comparison <- bundle$comparison
  if (is.null(comparison) || is.null(comparison$ollama) || is.null(comparison$openai)) return(NULL)
  provider_card <- function(provider_bundle, label, premium = FALSE) {
    entry <- provider_bundle$agents[[agent_id]] %or_else% list()
    provider_complete <- provider_bundle$source %in% c("ollama", "openai")
    summary <- if (provider_complete) {
      as.character(entry$summary %or_else% "The provider completed without a displayable agent summary.")
    } else {
      paste(
        label,
        "did not generate model-authored prose.",
        as.character(provider_bundle$reason %or_else% "The deterministic scientific summary remains available.")
      )
    }
    usage <- provider_bundle$usage %or_else% list()
    tokens <- suppressWarnings(as.numeric(usage$total_tokens %or_else% NA_real_))
    cost <- suppressWarnings(as.numeric(provider_bundle$estimated_cost_usd %or_else% NA_real_))
    shiny::div(
      class = paste("provider-comparison-card", if (premium) "provider-comparison-premium" else "provider-comparison-local"),
      shiny::div(
        class = "provider-comparison-heading",
        shiny::span(if (premium) "PREMIUM" else "LOCAL"),
        shiny::strong(label),
        shiny::tags$small(provider_bundle$model %or_else% "Computed summary")
      ),
      shiny::div(
        class = "provider-comparison-metrics",
        shiny::span(paste(interpretation_word_count(summary), "words")),
        if (is.finite(provider_bundle$elapsed_seconds %or_else% NA_real_)) shiny::span(paste0(round(provider_bundle$elapsed_seconds, 1), " s")) else NULL,
        if (is.finite(tokens)) shiny::span(paste(format(tokens, big.mark = ","), "tokens")) else NULL,
        if (is.finite(cost)) shiny::span(paste0("~$", format(cost, digits = 3, nsmall = 3), " USD")) else NULL
      ),
      shiny::p(summary),
      shiny::tags$small(
        class = "provider-grounding-note",
        if (provider_complete) {
          "Named findings and observed-result bullets remain deterministic and use the same result digest in both modes."
        } else {
          paste("Provider state:", interpretation_source_badge(provider_bundle), "· No fallback text is presented as provider output.")
        }
      )
    )
  }
  shiny::tags$details(
    class = "provider-comparison-shell",
    open = "open",
    shiny::tags$summary("Ollama vs OpenAI · controlled interpretation comparison"),
    shiny::div(
      class = "provider-comparison-grid",
      provider_card(comparison$ollama, "Ollama", premium = FALSE),
      provider_card(comparison$openai, "OpenAI", premium = TRUE)
    ),
    shiny::p(class = "provider-comparison-disclaimer", paste("Primary report narrative:", toupper(comparison$primary_provider %or_else% "computed"), "· Longer prose is not automatically more scientifically correct."))
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
      shiny::tagList(
        interpretation_source_ui(bundle),
        shiny::tags$p(
          class = "interpretation-empty",
          paste(
            "No significant result rows are available for this agent.",
            "This can be a valid biological outcome depending on the input data and thresholds."
          )
        )
      )
    )
  }

  evidence_items <- lapply(entry$observed_results %or_else% entry$evidence %or_else% character(), shiny::tags$li)
  key_findings <- lapply(entry$key_findings %or_else% character(), shiny::tags$li)

  observed_layer <- shiny::tags$details(
    class = "interpretation-details interpretation-observed-results",
    open = "open",
    shiny::tags$summary("Observed results · deterministic"),
    shiny::div(
      class = "interpretation-section",
      shiny::p("Computed from the complete saved result table; local AI cannot alter these observations."),
      shiny::tags$ul(evidence_items)
    )
  )

  if (is_interpretation_progress_bundle(bundle)) {
    return(shiny::tagList(
      interpretation_source_ui(bundle),
      observed_layer,
      shiny::div(
        class = "interpretation-loading-card",
        shiny::div(class = "interpretation-loading-orbit", "AI"),
        shiny::div(
          shiny::h4("Building the integrated interpretation"),
          shiny::p("Reviewing all result rows, checking grounding, and organizing one integrated biological account. No provisional computed narrative is shown while generation is active.")
        )
      )
    ))
  }

  shiny::tagList(
    interpretation_source_ui(bundle),
    provider_comparison_agent_ui(bundle, agent_id),
    observed_layer,
    shiny::div(
      class = "interpretation-section interpretation-biological-context",
      shiny::h4("Integrated biological interpretation"),
      shiny::p(entry$summary),
      if (length(key_findings)) shiny::tagList(shiny::h4("Key findings"), shiny::tags$ul(key_findings)) else NULL,
      shiny::h4("Biological context"),
      shiny::p(entry$biological_context %or_else% "Additional biological context is not available.")
    ),
    shiny::div(
      class = "interpretation-section",
      shiny::h4("Cancer relevance"),
      shiny::p(entry$cancer_relevance)
    )
  )
}

deep_narrative_ui <- function(value) {
  lines <- strsplit(as.character(value %or_else% ""), "\\r?\\n", perl = TRUE)[[1L]]
  lines <- lines[nzchar(trimws(lines))]
  suppress_section <- FALSE
  rendered <- lapply(lines, function(line) {
    text <- trimws(gsub("\\*\\*", "", line, perl = TRUE))
    known_sections <- c(
      "integrated biological interpretation", "evidence convergence and distinctions",
      "candidate regulatory and network model", "result-grounded research hypotheses",
      "testable hypotheses", "validation and next analyses", "recommended next analyses",
      "interpretation limits", "interpretive limits", "interpretive boundaries"
    )
    heading <- grepl("^\\s*\\*\\*.*\\*\\*\\s*$", line, perl = TRUE) || tolower(text) %in% known_sections
    if (heading) {
      suppress_section <<- tolower(text) %in% c(
        "result-grounded research hypotheses", "testable hypotheses",
        "validation and next analyses", "recommended next analyses",
        "interpretation limits", "interpretive limits", "interpretive boundaries"
      )
      if (suppress_section) return(NULL)
    } else if (suppress_section) {
      return(NULL)
    }
    if (heading) return(shiny::h4(text))
    if (grepl("^[0-9]+\\.\\s+", text, perl = TRUE)) return(shiny::div(class = "deep-narrative-item", text))
    if (grepl("^[*+-]\\s+", text, perl = TRUE)) return(shiny::div(class = "deep-narrative-subitem", sub("^[*+-]\\s+", "", text, perl = TRUE)))
    shiny::p(text)
  })
  shiny::tagList(Filter(Negate(is.null), rendered))
}

consolidated_hypothesis_ui <- function(synthesis) {
  convergences <- utils::head(as.character(synthesis$convergences %or_else% character()), 3L)
  hubs <- utils::head(as.character(synthesis$hub_candidates %or_else% character()), 3L)
  signals <- if (length(convergences)) paste(convergences, collapse = "; ") else "the recurring enriched programs"
  anchors <- if (length(hubs)) paste(hubs, collapse = ", ") else "the highest-ranked network hubs"
  shiny::div(
    class = "cross-agent-hypothesis",
    shiny::h4("One consolidated result-grounded hypothesis"),
    shiny::p(paste0(
      "Working hypothesis: the recurring signals (", signals,
      ") define a coordinated gene program in which ", anchors,
      " are candidate network anchors. Test this single model in an independent biological contrast using targeted perturbation or orthogonal protein/RNA measurements, with prespecified positive and negative controls."
    )),
    shiny::tags$small("This is a result-grounded hypothesis for validation, not a causal or clinical conclusion.")
  )
}

build_cross_agent_synthesis_ui <- function(bundle) {
  synthesis <- bundle$synthesis
  if (is_interpretation_progress_bundle(bundle)) {
    return(div(
      class = "cross-agent-content",
      interpretation_source_ui(bundle),
      div(class = "interpretation-loading-card", h4("Integrating evidence across agents"), p("The local model is completing groundedness checks and cross-agent synthesis."))
    ))
  }
  convergences <- lapply(synthesis$convergences %or_else% character(), tags$li)
  hubs <- lapply(synthesis$hub_candidates %or_else% character(), tags$li)
  bridge_available <- isTRUE(synthesis$bridge$available)

  div(
    class = "cross-agent-content",
    interpretation_source_ui(bundle),
    div(class = "cross-agent-summary", h3(synthesis$title %or_else% "Integrated interpretation"), p(synthesis$summary)),
    div(
      class = "cross-agent-narrative",
      h4(if (!is.null(synthesis$deep_narrative)) "Deep integrated interpretation" else "Integrated biological model"),
      div(class = "deep-narrative-text", deep_narrative_ui(synthesis$deep_narrative %or_else% synthesis$integrated_interpretation %or_else% synthesis$summary))
    ),
    div(
      class = "cross-agent-grid",
      div(h4("Convergent signals"), tags$ul(convergences)),
      div(h4("Candidate regulatory network"), p(synthesis$regulatory_network)),
      div(h4("Hub candidates"), tags$ul(hubs)),
      div(
        class = if (bridge_available) "drug-pathway-bridge bridge-ready" else "drug-pathway-bridge",
        h4("Drug ↔ pathway exchange"),
        p(synthesis$drug_pathway_context)
      ),
      div(h4("Novelty and literature context"), p(synthesis$novelty_context)),
      consolidated_hypothesis_ui(synthesis)
    )
  )
}

start_local_interpretation_job <- function(data_by_agent, settings) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("Background interpretation requires the processx package.", call. = FALSE)
  }

  worker_path <- normalizePath(
    file.path(project_root, "shiny-app", "run_interpretation_worker.R"),
    mustWork = TRUE
  )
  helper_path <- normalizePath(
    file.path(project_root, "shiny-app", "interpretation.R"),
    mustWork = TRUE
  )
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

  settings <- normalise_interpretation_settings(settings)
  started_at <- Sys.time()
  timeout_budget <- switch(
    settings$provider,
    openai = settings$openai_timeout_seconds * 2,
    compare = settings$timeout_seconds * 2 + settings$openai_timeout_seconds * 2,
    settings$timeout_seconds * 2
  )
  watchdog_grace_seconds <- min(10, max(2, timeout_budget * 0.02))

  list(
    process = process,
    input_path = input_path,
    output_path = output_path,
    stdout_path = stdout_path,
    stderr_path = stderr_path,
    settings = settings,
    started_at = started_at,
    # One bounded revision pass may be requested when the first structured
    # draft is valid but underdeveloped.
    deadline_at = started_at + timeout_budget + watchdog_grace_seconds
  )
}

local_interpretation_job_expired <- function(job, now = Sys.time()) {
  !is.null(job$deadline_at) && isTRUE(now >= job$deadline_at)
}

read_local_interpretation_job_bundle <- function(job) {
  if (is.null(job$output_path) || !file.exists(job$output_path)) return(NULL)
  bundle <- tryCatch(readRDS(job$output_path), error = function(error) NULL)
  if (is_terminal_interpretation_bundle(bundle)) bundle else NULL
}

cleanup_local_interpretation_job <- function(job, terminate = FALSE) {
  if (is.null(job)) return(invisible(NULL))
  if (isTRUE(terminate) && !is.null(job$process)) {
    try({
      if (isTRUE(job$process$is_alive())) {
        job$process$kill_tree()
        job$process$wait(timeout = 1000)
      }
    }, silent = TRUE)
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
  ollama_settings = NULL,
  run_context = NULL,
  interpretation_job_factory = start_local_interpretation_job,
  interpretation_cache_path = completed_interpretation_cache_path
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
    settings <- normalise_interpretation_settings(settings)

    run_is_complete <- if (is.null(analysis_complete)) TRUE else isTRUE(analysis_complete())
    csv_paths <- vapply(keys, function(key) result_files[[key]]$csv, character(1))
    information <- file.info(csv_paths)
    cache_key <- paste(
      keys,
      file.exists(csv_paths),
      information$size,
      as.numeric(information$mtime),
      settings$enabled,
      settings$provider,
      settings$host,
      settings$model,
      settings$timeout_seconds,
      settings$num_predict,
      settings$openai_model,
      settings$openai_reasoning_effort,
      settings$openai_timeout_seconds,
      settings$openai_max_output_tokens,
      settings$openai_data_consent,
      settings$openai_key_available,
      interpretation_contract_version,
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
          "Scientific workers are still running. The selected interpretation mode starts automatically when the analyses finish."
        ),
        cache_key
      )
      return()
    }

    persisted <- read_completed_interpretation(cache_key, interpretation_cache_path)
    if (!is.null(persisted)) {
      publish_interpretation(persisted, cache_key)
      return()
    }

    if (!isTRUE(settings$enabled)) {
      publish_interpretation(
        build_rule_interpretation_bundle(
          exchanges,
          "AI interpretation is disabled; the deterministic full-result summary is shown."
        ),
        cache_key
      )
      return()
    }

    validation_error <- validate_interpretation_settings(settings)
    if (!is.null(validation_error)) {
      bundle <- build_provider_failure_bundle(
        exchanges,
        simpleError(validation_error),
        settings,
        reason = paste("The selected interpretation mode was not started:", validation_error)
      )
      persist_completed_interpretation(cache_key, bundle, interpretation_cache_path)
      publish_interpretation(bundle, cache_key)
      return()
    }

    loading <- build_provider_progress_bundle(
      exchanges,
      state = "loading",
      settings = settings
    )
    publish_interpretation(loading, cache_key)

    interpretation_job$value <- tryCatch(
      {
        job <- interpretation_job_factory(data_by_agent, settings)
        job$key <- cache_key
        job
      },
      error = function(error) {
        fallback <- build_provider_failure_bundle(exchanges, error, settings)
        persist_completed_interpretation(cache_key, fallback, interpretation_cache_path)
        publish_interpretation(fallback, cache_key)
        NULL
      }
    )
  })

  finalize_interpretation_job <- function() {
    job <- interpretation_job$value
    if (is.null(job)) return(invisible(FALSE))

    is_alive <- isTRUE(tryCatch(job$process$is_alive(), error = function(error) FALSE))
    bundle <- read_local_interpretation_job_bundle(job)
    if (!is.null(bundle)) {
      completed_key <- job$key
      cleanup_local_interpretation_job(job, terminate = is_alive)
      interpretation_job$value <- NULL
      if (identical(completed_key, interpretation_cache$key)) {
        persist_completed_interpretation(completed_key, bundle, interpretation_cache_path)
        publish_interpretation(bundle, completed_key)
      }
      return(invisible(TRUE))
    }

    if (is_alive) {
      current <- interpretation_cache$value
      if (local_interpretation_job_expired(job)) {
        timeout_error <- simpleError("Interpretation timed out at the configured time limit.")
        bundle <- build_provider_failure_bundle(
          current$exchanges %or_else% list(),
          timeout_error,
          job$settings,
          source = if (identical(job$settings$provider, "openai")) "openai_timeout" else if (identical(job$settings$provider, "compare")) "comparison_timeout" else "ollama_timeout"
        )
        completed_key <- job$key
        cleanup_local_interpretation_job(job, terminate = TRUE)
        interpretation_job$value <- NULL
        if (identical(completed_key, interpretation_cache$key)) {
          persist_completed_interpretation(completed_key, bundle, interpretation_cache_path)
          publish_interpretation(bundle, completed_key)
        }
        return(invisible(TRUE))
      }
      if (
        identical(job$key, interpretation_cache$key) &&
        identical(current$source, "loading")
      ) {
        generating <- build_provider_progress_bundle(
          current$exchanges %or_else% list(),
          state = "generating",
          settings = job$settings
        )
        publish_interpretation(generating, job$key)
      }
      return(invisible(FALSE))
    }

    if (is.null(bundle)) {
      current <- interpretation_cache$value
      exchanges <- current$exchanges %or_else% list()
      bundle <- build_provider_failure_bundle(
        exchanges,
        simpleError("The background local-model process ended without a valid response."),
        job$settings,
        source = if (identical(job$settings$provider, "openai")) "openai_error" else if (identical(job$settings$provider, "compare")) "comparison_error" else "ollama_error",
        reason = "The background interpretation process ended without a valid response. Scientific results remain available with the computed scientific summary."
      )
    }

    completed_key <- job$key
    cleanup_local_interpretation_job(job)
    interpretation_job$value <- NULL
    if (identical(completed_key, interpretation_cache$key)) {
      persist_completed_interpretation(completed_key, bundle, interpretation_cache_path)
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
      context <- if (is.null(run_context)) list() else run_context()
      build_combined_html_report(file, current, active_keys(), run_context = context)
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
      context <- if (is.null(run_context)) list() else run_context()
      create_results_bundle(file, current, active_keys(), run_context = context)
    },
    contentType = "application/zip"
  )

  for (agent_key in names(agent_titles)) {
    local({
      key <- agent_key
      title <- agent_titles[[key]]
      csv_path <- result_files[[key]]$csv
      plot_path <- result_files[[key]]$plot

      visual_summary <- reactive({
        result_signature()
        result_visual_summary(safe_result_csv(csv_path), key)
      })

      visual_interactions <- reactive({
        result_signature()
        if (!identical(key, "string")) return(NULL)
        safe_result_csv(result_files$string$interactions)
      })

      output[[paste0(key, "_interactive_bar")]] <- plotly::renderPlotly({
        professional_bar_plot(visual_summary(), paste(title, "evidence profile"))
      })

      if (identical(key, "string")) {
        output[[paste0(key, "_interactive_3d")]] <- plotly::renderPlotly({
          evidence_3d_plot(visual_summary(), title, visual_interactions())
        })

        output[[paste0(key, "_network_guide")]] <- renderUI({
          string_network_guide_ui(visual_summary(), visual_interactions())
        })
      }

      output[[paste0(key, "_interactive_bar_modal")]] <- plotly::renderPlotly({
        professional_bar_plot(visual_summary(), paste(title, "evidence profile"))
      })

      observeEvent(input[[paste0("open_plot_", key)]], {
        if (is.null(visual_summary())) {
          showNotification(paste(title, "visualization is not available yet."), type = "warning")
          return()
        }

        showModal(
          modalDialog(
            title = paste(title, "interactive evidence profile"),
            plotly::plotlyOutput(paste0(key, "_interactive_bar_modal"), height = "70vh"),
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
        summary <- visual_summary()
        if (is.null(summary) || !requireNamespace("plotly", quietly = TRUE)) {
          result_data <- safe_result_csv(csv_path)
          if (is_existing_result_file(plot_path)) {
            plot_source <- image_data_uri(plot_path)
            if (!is.null(plot_source)) return(tags$img(src = plot_source, alt = paste(title, "analysis visualization"), class = "result-plot-image"))
          }
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

        bar_view <- plotly::plotlyOutput(paste0(key, "_interactive_bar"), height = "470px")
        if (!identical(key, "string")) return(bar_view)

        network <- evidence_3d_plot(summary, title, visual_interactions())
        if (is.null(network)) {
          return(tagList(
            bar_view,
            div(
              class = "result-empty-state string-network-unavailable",
              h4("Connected STRING network unavailable"),
              p("No retrieved interaction edges connected the displayed hub proteins. The bar view and full CSV remain available; no artificial rank-order lines were drawn.")
            )
          ))
        }

        tabsetPanel(
          type = "pills",
          tabPanel("Professional bar view", bar_view),
          tabPanel(
            "3D interaction network",
            div(
              class = "string-app-network-layout",
              div(
                class = "string-app-network-chart",
                p(
                  class = "visualization-method-note",
                  "Drag to rotate and scroll to zoom. Balls are hub proteins; lines are retrieved STRING associations."
                ),
                plotly::plotlyOutput(paste0(key, "_interactive_3d"), height = "540px")
              ),
              uiOutput(paste0(key, "_network_guide"))
            ),
          )
        )
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
