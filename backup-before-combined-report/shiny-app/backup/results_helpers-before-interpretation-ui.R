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
          uiOutput(paste0(key, "_result_status"))
        ),

        uiOutput(paste0(key, "_result_plot"))
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
        class = "results-ready-badge",
        span(class = "results-ready-dot"),
        "Saved outputs"
      )
    ),

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

register_results_server <- function(output, session) {

  agent_titles <- c(
    go = "GO",
    kegg = "KEGG",
    gsva = "GSVA",
    chea = "ChEA"
  )

  for (agent_key in names(agent_titles)) {
    local({
      key <- agent_key
      title <- agent_titles[[key]]
      csv_path <- result_files[[key]]$csv
      plot_path <- result_files[[key]]$plot

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
