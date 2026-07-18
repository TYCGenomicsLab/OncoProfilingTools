project_root <- if (basename(getwd()) == "shiny-app") {
  normalizePath("..", mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

analysis_output_root <- file.path(project_root, "output")

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
  )
)

if (dir.exists(analysis_output_root)) {
  shiny::addResourcePath(
    prefix = "analysis-output",
    directoryPath = analysis_output_root
  )
}

safe_result_csv <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble())
  }

  tryCatch(
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
}


html_escape_value <- function(value) {
  htmltools::htmlEscape(as.character(if (is.null(value)) "" else value))
}

image_data_uri <- function(path) {
  if (!file.exists(path) || !requireNamespace("base64enc", quietly = TRUE)) {
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

result_summary_rows <- function() {
  keys <- names(result_files)

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
        Plot = if (file.exists(plot_path)) "Available" else "Not available",
        stringsAsFactors = FALSE
      )
    })
  )
}

build_combined_html_report <- function(destination) {
  titles <- c(go = "GO", kegg = "KEGG", gsva = "GSVA", chea = "ChEA")
  generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  summary_data <- result_summary_rows()

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

  sections <- vapply(names(result_files), function(key) {
    csv_path <- result_files[[key]]$csv
    plot_path <- result_files[[key]]$plot
    data <- safe_result_csv(csv_path)
    title <- titles[[key]]
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
      preview <- utils::head(data, 25)
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
      "<p class='section-copy'>", html_escape_value(
        switch(key,
          go = "Gene Ontology biological-process enrichment.",
          kegg = "KEGG pathway enrichment.",
          gsva = "Hallmark pathway activity across samples.",
          chea = "Transcription-factor enrichment from ChEA."
        )
      ), "</p>",
      image_html,
      "<h3>Top results</h3>", table_html,
      "</section>"
    )
  }, character(1))

  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>OncoProfiling Combined Analysis Report</title>",
    "<style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;background:#07111f;color:#e8eef7;line-height:1.55}",
    ".page{max-width:1180px;margin:auto;padding:40px 28px 80px}",
    ".hero,section{background:#0d1b2d;border:1px solid #20324a;border-radius:18px;padding:28px;margin-bottom:24px;box-shadow:0 18px 50px rgba(0,0,0,.24)}",
    "h1{font-size:34px;margin:0 0 8px}h2{font-size:25px;margin-top:0}h3{margin-top:28px}",
    ".muted,.note,.section-copy{color:#a9b8ca}.badge{display:inline-block;padding:6px 10px;border-radius:999px;background:#123f38;color:#8df1cf;font-size:13px;font-weight:700}",
    "figure{margin:22px 0}img{display:block;max-width:100%;max-height:720px;margin:auto;border-radius:12px;background:white}figcaption{text-align:center;color:#a9b8ca;margin-top:8px}",
    ".table-wrap{overflow:auto;border:1px solid #263a54;border-radius:12px}table{border-collapse:collapse;width:100%;font-size:13px;background:#0a1727}th,td{padding:10px 12px;border-bottom:1px solid #1d3048;text-align:left;vertical-align:top;white-space:nowrap}th{background:#13243a;position:sticky;top:0}.empty{padding:20px;border:1px dashed #38506e;border-radius:12px;color:#a9b8ca}",
    "</style></head><body><main class='page'>",
    "<div class='hero'><span class='badge'>COMBINED REPORT</span>",
    "<h1>OncoProfiling Multi-Agent Analysis</h1>",
    "<p class='muted'>Generated ", html_escape_value(generated_at), ". This report combines GO, KEGG, GSVA, and ChEA outputs into one portable HTML document.</p>",
    "<h3>Run summary</h3><div class='table-wrap'><table><thead><tr><th>Agent</th><th>Status</th><th>Rows</th><th>Plot</th></tr></thead><tbody>", summary_rows, "</tbody></table></div></div>",
    paste0(sections, collapse = "\n"),
    "</main></body></html>"
  )

  writeLines(html, destination, useBytes = TRUE)
}

create_results_bundle <- function(destination) {
  temp_dir <- tempfile("oncoprofiling-results-")
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  report_path <- file.path(temp_dir, "OncoProfiling_Combined_Report.html")
  build_combined_html_report(report_path)

  copied <- c(report_path)
  for (key in names(result_files)) {
    for (kind in c("csv", "plot")) {
      source_path <- result_files[[key]][[kind]]
      if (file.exists(source_path)) {
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

        uiOutput(paste0(key, "_result_plot")),

        div(
          class = "biological-interpretation-card",

          div(
            class = "interpretation-card-heading",
            div(class = "interpretation-card-icon", "AI"),
            div(
              span(class = "interpretation-card-label", "BIOLOGICAL INTERPRETATION"),
              h4(paste(title, "Agent Summary"))
            )
          ),

          uiOutput(paste0(key, "_interpretation"))
        )
      ),

      div(
        class = "result-table-panel result-inner-card",

        div(
          class = "result-panel-heading",
          div(
            h3(paste(title, "Results")),
            p("Search, filter, and inspect the generated result records.")
          ),

          downloadButton(
            outputId = paste0("download_", key),
            label = "Download CSV",
            class = "result-download-button"
          )
        ),

        DTOutput(paste0(key, "_result_table"))
      )
    )
  )
}

results_center_ui <- function() {
  div(
    id = "results-center",
    class = "results-center glass-card",

    div(
      class = "results-center-heading",

      div(
        div(class = "section-label", "STEP 3"),
        h2("Analysis Results Center"),
        p(
          class = "section-description",
          paste(
            "Explore generated outputs from GO, KEGG, GSVA, and ChEA.",
            "Each agent has an independent visualization and result table."
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

    tabsetPanel(
      id = "results_tabs",
      type = "pills",

      result_tab_ui(
        "go",
        "GO",
        "Enriched Gene Ontology biological-process terms."
      ),

      result_tab_ui(
        "kegg",
        "KEGG",
        "Significantly enriched molecular pathways."
      ),

      result_tab_ui(
        "gsva",
        "GSVA",
        "Hallmark pathway activity across biological samples."
      ),

      result_tab_ui(
        "chea",
        "ChEA",
        "Candidate transcription-factor regulators."
      )
    )
  )
}


build_agent_interpretation <- function(agent_id, data) {

  if (is.null(data) || nrow(data) == 0) {
    return(
      tags$p(
        class = "interpretation-empty",
        paste(
          "No significant result rows are available for this agent.",
          "This can be a valid biological outcome depending on",
          "the input data and significance thresholds."
        )
      )
    )
  }

  find_first_column <- function(candidates) {
    matched <- candidates[candidates %in% names(data)]

    if (length(matched) == 0) {
      return(NULL)
    }

    matched[[1]]
  }

  description_column <- find_first_column(
    c(
      "Description",
      "description",
      "Term",
      "term",
      "Pathway",
      "pathway",
      "GeneSet",
      "gene_set"
    )
  )

  adjusted_p_column <- find_first_column(
    c(
      "p.adjust",
      "Adjusted.P.value",
      "Adjusted.P.value.",
      "adj_p_value",
      "FDR",
      "qvalue"
    )
  )

  score_column <- find_first_column(
    c(
      "Combined.Score",
      "Combined Score",
      "NES",
      "enrichmentScore",
      "FoldEnrichment",
      "RichFactor",
      "Count"
    )
  )

  top_term <- if (!is.null(description_column)) {
    as.character(data[[description_column]][[1]])
  } else {
    "the highest-ranked result"
  }

  adjusted_p_text <- ""

  if (!is.null(adjusted_p_column)) {
    value <- suppressWarnings(
      as.numeric(data[[adjusted_p_column]][[1]])
    )

    if (!is.na(value)) {
      adjusted_p_text <- paste0(
        " Its adjusted significance value is ",
        format(
          value,
          scientific = value < 0.001,
          digits = 3
        ),
        "."
      )
    }
  }

  score_text <- ""

  if (!is.null(score_column)) {
    value <- suppressWarnings(
      as.numeric(data[[score_column]][[1]])
    )

    if (!is.na(value)) {
      score_text <- paste0(
        " The leading ranking metric is ",
        format(
          round(value, 3),
          big.mark = ",",
          trim = TRUE
        ),
        "."
      )
    }
  }

  interpretation <- switch(
    agent_id,

    go = paste0(
      "The strongest Gene Ontology biological-process signal is ",
      top_term,
      ". This suggests that the uploaded gene profile is most strongly ",
      "associated with this coordinated biological program.",
      adjusted_p_text,
      score_text
    ),

    kegg = paste0(
      "The highest-ranked KEGG pathway is ",
      top_term,
      ". This indicates that the input genes may converge on this ",
      "molecular pathway or disease-associated signaling system.",
      adjusted_p_text,
      score_text
    ),

    gsva = paste0(
      "The leading GSVA Hallmark signal is ",
      top_term,
      ". This represents elevated or distinctive pathway activity ",
      "across the analyzed samples rather than simple gene-list overlap.",
      adjusted_p_text,
      score_text
    ),

    chea = paste0(
      "The leading ChEA regulatory signal is ",
      top_term,
      ". This suggests that the corresponding transcription factor ",
      "may regulate a meaningful portion of the submitted gene program.",
      adjusted_p_text,
      score_text
    ),

    paste0(
      "The highest-ranked result is ",
      top_term,
      ".",
      adjusted_p_text,
      score_text
    )
  )

  tags$p(
    HTML(
      paste0(
        "<span class='interpretation-highlight'>",
        htmltools::htmlEscape(top_term),
        "</span>. ",
        htmltools::htmlEscape(
          sub(
            paste0("^.*?\\Q", top_term, "\\E\\. "),
            "",
            interpretation,
            perl = TRUE
          )
        )
      )
    )
  )
}

register_results_server <- function(input, output, session) {

  agent_titles <- c(
    go = "GO",
    kegg = "KEGG",
    gsva = "GSVA",
    chea = "ChEA"
  )

  output$analysis_completion_banner <- renderUI({
    summary_data <- result_summary_rows()
    completed <- sum(summary_data$Status != "Not generated")
    total_rows <- sum(summary_data$Rows)

    div(
      class = paste(
        "analysis-completion-banner",
        if (completed == 4) "analysis-completion-complete" else "analysis-completion-partial"
      ),
      div(
        class = "analysis-completion-icon",
        if (completed == 4) "✓" else "↻"
      ),
      div(
        h4(if (completed == 4) "Four-agent analysis outputs are ready" else "Analysis outputs are partially available"),
        p(
          paste0(
            completed, " of 4 agents generated output files · ",
            format(total_rows, big.mark = ","), " total result rows"
          )
        )
      )
    )
  })

  output$download_combined_report <- downloadHandler(
    filename = function() {
      paste0("OncoProfiling_Combined_Report_", format(Sys.Date(), "%Y%m%d"), ".html")
    },
    content = function(file) {
      build_combined_html_report(file)
    },
    contentType = "text/html"
  )

  output$download_all_results <- downloadHandler(
    filename = function() {
      paste0("OncoProfiling_All_Results_", format(Sys.Date(), "%Y%m%d"), ".zip")
    },
    content = function(file) {
      create_results_bundle(file)
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
        if (!file.exists(plot_path)) {
          showNotification(paste(title, "visualization is not available yet."), type = "warning")
          return()
        }

        relative_plot_path <- sub(
          paste0("^", analysis_output_root, "/?"),
          "",
          plot_path
        )

        showModal(
          modalDialog(
            title = paste(title, "Visualization"),
            tags$img(
              src = paste0(
                "analysis-output/", relative_plot_path,
                "?v=", as.numeric(file.info(plot_path)$mtime)
              ),
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
        result_data <- safe_result_csv(csv_path)

        build_agent_interpretation(
          agent_id = key,
          data = result_data
        )
      })

      output[[paste0(key, "_result_plot")]] <- renderUI({
        if (!file.exists(plot_path)) {
          message_text <- if (
            key == "kegg" &&
            file.exists(csv_path) &&
            nrow(safe_result_csv(csv_path)) == 0
          ) {
            paste(
              "KEGG completed successfully, but no pathways",
              "passed the current significance threshold."
            )
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

        tags$img(
          src = paste0(
            "analysis-output/",
            relative_plot_path,
            "?v=",
            as.numeric(file.info(plot_path)$mtime)
          ),
          alt = paste(title, "analysis visualization"),
          class = "result-plot-image"
        )
      })

      output[[paste0(key, "_result_table")]] <- renderDT({
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

        DT::datatable(
          result_data,
          rownames = FALSE,
          filter = "top",
          options = list(
            pageLength = 10,
            lengthMenu = c(5, 10, 25, 50),
            scrollX = TRUE,
            autoWidth = TRUE,
            dom = "tip"
          ),
          class = "stripe hover compact"
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
    })
  }
}
