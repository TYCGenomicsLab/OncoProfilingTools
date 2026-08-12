options(shiny.maxRequestSize = 1024 * 1024^2)

suppressPackageStartupMessages({
  library(shiny)
  library(readr)
  library(DT)
})

source("real_pipeline.R", local = TRUE)
source("run_real_agents.R", local = TRUE)
source("interpretation.R", local = TRUE)
source("results_helpers.R", local = TRUE)

# Semantic HTML helpers are not exported as top-level Shiny functions.
header <- tags$header
main <- tags$main
section <- tags$section
footer <- tags$footer

module_meta <- list(
  go = list(title = "GO", subtitle = "Biological process enrichment", input = "Gene list"),
  kegg = list(title = "KEGG", subtitle = "Molecular pathway enrichment", input = "Gene list"),
  reactome = list(title = "Reactome", subtitle = "Curated pathway enrichment", input = "Gene list"),
  wikipathways = list(title = "WikiPathways", subtitle = "Community pathway enrichment", input = "Gene list"),
  string = list(title = "STRING", subtitle = "Protein interaction hubs", input = "Gene list"),
  hallmark = list(title = "Hallmark", subtitle = "Cancer hallmark enrichment", input = "Gene list"),
  chea = list(title = "ChEA", subtitle = "Transcription factor enrichment", input = "Gene list"),
  gsva = list(title = "GSVA", subtitle = "Sample-level pathway activity", input = "Expression matrix"),
  immune = list(title = "Immune Deconvolution", subtitle = "Immune cell composition", input = "Expression matrix"),
  drug = list(title = "Drug Sensitivity", subtitle = "Ranked compound responses", input = "Drug-response table")
)

read_dataset <- function(path, name) {
  extension <- tolower(tools::file_ext(name))
  switch(
    extension,
    csv = readr::read_csv(path, show_col_types = FALSE, name_repair = "unique"),
    tsv = readr::read_tsv(path, show_col_types = FALSE, name_repair = "unique"),
    txt = readr::read_tsv(path, show_col_types = FALSE, name_repair = "unique"),
    stop("Upload a CSV, TSV, or TXT file.", call. = FALSE)
  )
}

dataset_profile <- function(data) {
  empty <- list(
    gene = FALSE, expression = FALSE, drug = FALSE,
    genes = character(), analysis_genes = character(), gene_column = NULL,
    selection_note = NULL, selection_error = NULL
  )
  if (is.null(data) || nrow(data) == 0L || ncol(data) == 0L) return(empty)

  normalized <- normalise_column_name(names(data))
  gene_candidates <- c("genesymbol", "hugosymbol", "gene", "genes", "symbol", "genename", "hgncsymbol")
  gene_positions <- match(gene_candidates, normalized, nomatch = 0L)
  gene_positions <- gene_positions[gene_positions > 0L]
  gene_column <- if (length(gene_positions)) names(data)[gene_positions[[1]]] else NULL
  genes <- if (!is.null(gene_column)) clean_gene_symbols(data[[gene_column]]) else character()

  numeric_count <- sum(vapply(data, function(x) {
    values <- suppressWarnings(as.numeric(as.character(x)))
    mean(!is.na(values)) >= 0.8
  }, logical(1)))
  # DEG summaries such as gene_expression_CMS4.csv have a few numeric
  # statistics but are not sample-by-gene matrices. Require a wider matrix
  # before offering GSVA or immune deconvolution.
  expression <- ncol(data) >= 10L && length(genes) >= 10L && numeric_count >= 2L

  compound <- any(normalized %in% c("compound", "compoundname", "drug", "drugname", "treatment"))
  response <- any(normalized %in% c("ic50", "auc", "viability", "sensitivity", "response", "lnic50"))
  wide_prism <- any(grepl("BRD:", names(data), fixed = TRUE))

  selection <- NULL
  selection_error <- NULL
  if (length(genes) >= 2L) {
    selection <- tryCatch(
      prepare_gene_input(data),
      error = function(error) {
        selection_error <<- conditionMessage(error)
        NULL
      }
    )
  }
  analysis_genes <- if (is.null(selection)) character() else selection$genes

  list(
    gene = length(analysis_genes) >= 2L,
    expression = expression,
    drug = (compound && response) || wide_prism,
    genes = genes,
    analysis_genes = analysis_genes,
    gene_column = gene_column,
    selection_note = if (is.null(selection)) NULL else selection$selection_note,
    selection_error = selection_error
  )
}

module_card <- function(key, meta) {
  div(
    id = paste0("module-card-", key),
    class = "module-card",
    checkboxInput(paste0("module_", key), label = NULL, value = FALSE),
    div(
      class = "module-copy",
      div(class = "module-title-row", strong(meta$title), span(id = paste0("compat-", key), class = "compat-badge", "Awaiting data")),
      span(class = "module-subtitle", meta$subtitle),
      span(class = "module-input", meta$input)
    )
  )
}

module_selector_group <- function(group_id, group) {
  keys <- intersect(group$keys, names(module_meta))
  div(
    class = paste("module-selector-group", paste0("module-selector-", group_id)),
    div(
      class = "module-group-heading",
      div(h3(group$title), p(group$description)),
      span(paste(length(keys), if (length(keys) == 1L) "agent" else "agents"))
    ),
    div(class = "module-grid", Map(module_card, keys, module_meta[keys]))
  )
}

progress_row <- function(key, meta) {
  div(
    id = paste0("progress-", key), class = "progress-row progress-idle",
    span(class = "progress-indicator", "·"),
    div(strong(meta$title), span(class = "progress-message", "Waiting"))
  )
}

ui <- fluidPage(
  tags$head(
    tags$title("OncoProfilingTools"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$link(rel = "stylesheet", href = "styles.css?v=release-candidate-4"),
    tags$script(src = "status.js?v=results-center-polish-2")
  ),
  div(
    class = "app-shell",
    header(
      class = "site-header",
      div(class = "brand-mark", "OP"),
      div(class = "brand-copy", h1("OncoProfilingTools"), p("Agentic AI Platform for Cancer Pharmacogenomics")),
      div(class = "system-status", span(), "Local-only workspace")
    ),

    main(
      class = "dashboard",
      section(
        class = "hero-panel",
        div(
          class = "hero-copy",
          span(class = "eyebrow", "PRIVATE BIOMARKER DISCOVERY WORKSPACE"),
          h2("Analyze once. Explore every biological signal."),
          p("Validated enrichment, network, immune, and pharmacogenomic workflows with local interpretation and exportable evidence.")
        ),
        div(class = "hero-facts", div(strong(as.character(length(module_meta))), span("Analysis modules")), div(strong("Real"), span("Scientific backends")))
      ),

      div(
        class = "workflow-grid",
        section(
          class = "panel input-panel",
          div(class = "section-heading", span(class = "step", "01"), div(h2("Upload dataset"), p("Gene list, expression matrix, or drug-response table"))),
          div(
            class = "upload-zone",
            div(class = "upload-icon", "CSV"),
            div(h3("Choose analysis file"), p("CSV, TSV, or TXT · up to 1 GB")),
            fileInput("dataset", label = NULL, accept = c(".csv", ".tsv", ".txt"), buttonLabel = "Choose file", placeholder = "No file selected")
          ),
          uiOutput("dataset_status"),
          uiOutput("dataset_metrics"),
          uiOutput("gene_selection_note")
        ),

        section(
          class = "panel module-panel",
          div(class = "section-heading", span(class = "step", "02"), div(h2("Analysis modules"), p("Only scientifically compatible modules will execute"))),
          lapply(names(agent_groups), function(group_id) module_selector_group(group_id, agent_groups[[group_id]])),
          tags$details(
            class = "ollama-settings",
            tags$summary("Local biological interpretation settings"),
            div(
              class = "ollama-settings-body",
              checkboxInput("ollama_enabled", "Use local Ollama when available", value = TRUE),
              div(
                class = "ollama-settings-grid",
                textInput("ollama_host", "Ollama host", value = "http://127.0.0.1:11434"),
                textInput("ollama_model", "Model", value = "llama3.1:8b"),
                numericInput("ollama_timeout", "Generation timeout (seconds)", value = 180, min = 30, max = 600, step = 30)
              ),
              p("All prompts stay on this computer. Only localhost, 127.0.0.1, or [::1] are accepted. Interpretation runs in a background process and safely falls back if the local model cannot respond.")
            )
          ),
          actionButton("run_analysis", "Run analysis", class = "primary-button", icon = icon("play")),
          p(class = "run-note", "Modules requiring a different input type are retained in the summary as not executed.")
        )
      ),

      section(
        id = "progress-section", class = "panel progress-panel",
        div(class = "section-heading", span(class = "step", "03"), div(h2("Analysis progress"), p("Independent background workers report live status"))),
        div(class = "progress-layout", div(class = "progress-list", Map(progress_row, names(module_meta), module_meta)), uiOutput("run_overview"))
      ),

      uiOutput("final_summary"),
      uiOutput("results_center")
    ),
    footer(class = "site-footer", span("OncoProfilingTools · Research use only"), span("Biological findings require expert validation"))
  )
)

server <- function(input, output, session) {
  data_value <- reactiveVal(NULL)
  data_error <- reactiveVal(NULL)
  workers_value <- reactiveVal(list())
  results_value <- reactiveVal(list())
  selected_value <- reactiveVal(character())
  run_started <- reactiveVal(NULL)
  run_finished <- reactiveVal(NULL)
  running <- reactiveVal(FALSE)

  ollama_settings <- reactive({
    list(
      enabled = isTRUE(input$ollama_enabled),
      host = input$ollama_host %or_else% "http://127.0.0.1:11434",
      model = input$ollama_model %or_else% "llama3.1:8b",
      timeout_seconds = input$ollama_timeout %or_else% 180,
      num_predict = 1800
    )
  })

  analysis_complete <- reactive({
    !running() && !is.null(run_finished())
  })

  stop_active_workers <- function() {
    workers <- shiny::isolate(workers_value())
    if (!length(workers)) return(invisible(NULL))

    for (worker in workers) {
      process <- worker$process
      if (is.null(process)) next
      try({
        if (isTRUE(process$is_alive())) process$kill()
      }, silent = TRUE)
    }

    invisible(NULL)
  }

  clear_previous_artifacts <- function() {
    # Remove every result file known to the Results Center. This prevents
    # plots, tables and downloadable reports from a previous dataset/run
    # from appearing in the current analysis.
    paths <- c(
      unique(unlist(lapply(result_files, function(x) unname(unlist(x))), use.names = FALSE)),
      completed_interpretation_cache_path
    )
    paths <- paths[!is.na(paths) & nzchar(paths)]
    existing <- paths[file.exists(paths)]
    if (length(existing)) unlink(existing, recursive = TRUE, force = TRUE)

    # Remove background-worker inputs/results/logs from earlier runs.
    if (dir.exists(real_pipeline_runtime_dir)) {
      stale <- list.files(real_pipeline_runtime_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
      if (length(stale)) unlink(stale, recursive = TRUE, force = TRUE)
    }

    invisible(NULL)
  }

  reset_analysis_state <- function(clear_files = TRUE) {
    stop_active_workers()
    running(FALSE)
    workers_value(list())
    results_value(list())
    selected_value(character())
    run_started(NULL)
    run_finished(NULL)
    if (isTRUE(clear_files)) clear_previous_artifacts()
    session$sendCustomMessage("reset-analysis-ui", list())
    invisible(NULL)
  }

  # Never expose output files left behind by an earlier app session.
  clear_previous_artifacts()

  session$onSessionEnded(function() {
    stop_active_workers()
  })

  output$results_center <- renderUI({
    if (is.null(run_started())) return(NULL)
    results_center_ui(selected_value())
  })

  results_controller <- register_results_server(
    input,
    output,
    session,
    active_agents = selected_value,
    analysis_complete = analysis_complete,
    ollama_settings = ollama_settings
  )

  profile <- reactive(dataset_profile(data_value()))

  compatibility <- reactive({
    p <- profile()
    c(
      go = p$gene, kegg = p$gene, reactome = p$gene, wikipathways = p$gene, string = p$gene, hallmark = p$gene, chea = p$gene,
      gsva = p$expression, immune = p$expression, drug = p$drug
    )
  })

  observeEvent(input$dataset, {
    reset_analysis_state(clear_files = TRUE)
    data_error(NULL)
    data_value(NULL)
    tryCatch({
      loaded <- read_dataset(input$dataset$datapath, input$dataset$name)
      if (!nrow(loaded) || !ncol(loaded)) stop("The uploaded file is empty.")
      data_value(loaded)
    }, error = function(error) data_error(conditionMessage(error)))
  }, ignoreInit = TRUE)

  observe({
    compatible <- compatibility()
    for (key in names(module_meta)) {
      session$sendCustomMessage("module-compatibility", list(
        key = key,
        compatible = isTRUE(compatible[[key]]),
        ready = !is.null(data_value())
      ))
    }
  })

  output$dataset_status <- renderUI({
    if (!is.null(data_error())) return(div(class = "dataset-alert dataset-error", data_error()))
    if (is.null(input$dataset)) return(div(class = "dataset-alert", "Select a file to validate its analysis compatibility."))
    if (is.null(data_value())) return(div(class = "dataset-alert", "Reading dataset…"))
    div(class = "dataset-alert dataset-ready", span("✓"), strong(input$dataset$name), " validated")
  })

  output$dataset_metrics <- renderUI({
    data <- data_value()
    if (is.null(data)) return(NULL)
    p <- profile()
    div(
      class = "dataset-metrics",
      div(span("Rows"), strong(format(nrow(data), big.mark = ","))),
      div(span("Columns"), strong(format(ncol(data), big.mark = ","))),
      div(
        span("Analysis genes"),
        strong(paste0(format(length(p$analysis_genes), big.mark = ","), " / ", format(length(p$genes), big.mark = ",")))
      ),
      div(span("Compatible"), strong(sum(compatibility())))
    )
  })

  output$gene_selection_note <- renderUI({
    if (is.null(data_value())) return(NULL)
    p <- profile()
    if (!is.null(p$selection_error)) {
      return(div(class = "gene-selection-note gene-selection-warning", strong("Gene-list input needs attention"), p$selection_error))
    }
    if (is.null(p$selection_note)) return(NULL)
    div(
      class = "gene-selection-note",
      strong("Automatic gene-set preparation"),
      p$selection_note,
      tags$small("For directional biology, upload up- and down-regulated signatures as separate runs.")
    )
  })

  send_progress <- function(key, status, message) {
    session$sendCustomMessage("analysis-progress", list(key = key, status = status, message = message))
  }

  observeEvent(input$run_analysis, {
    data <- data_value()
    if (is.null(data)) {
      showNotification("Upload a valid dataset first.", type = "warning")
      return()
    }
    if (running()) return()

    requested <- vapply(names(module_meta), function(key) isTRUE(input[[paste0("module_", key)]]), logical(1))
    selected <- names(requested)[requested & compatibility()]
    if (!length(selected)) {
      showNotification("Select at least one compatible module.", type = "warning")
      return()
    }

    # A new run always starts from an empty state and empty output folders.
    # This guarantees that downloadable reports contain only this run.
    reset_analysis_state(clear_files = TRUE)
    selected_value(selected)
    run_started(Sys.time())
    running(TRUE)

    for (key in names(module_meta)) {
      if (key %in% selected) send_progress(key, "queued", "Queued") else send_progress(key, "skipped", "Not executed")
    }

    analysis_run <- tryCatch(create_real_analysis_run(data), error = function(error) {
      running(FALSE)
      showNotification(conditionMessage(error), type = "error")
      NULL
    })
    if (is.null(analysis_run)) return()

    workers <- list()
    for (key in selected) {
      send_progress(key, "running", paste("Running", module_meta[[key]]$title, "…"))
      workers[[key]] <- tryCatch(
        start_real_agent_worker(key, analysis_run, 0.05),
        error = function(error) list(agent = key, launch_error = conditionMessage(error), handled = FALSE)
      )
      workers[[key]]$handled <- FALSE
    }
    workers_value(workers)
  })

  observe({
    if (!running()) return()
    invalidateLater(600, session)
    workers <- workers_value()
    results <- results_value()
    changed <- FALSE

    for (key in names(workers)) {
      worker <- workers[[key]]
      if (isTRUE(worker$handled)) next
      if (!is.null(worker$launch_error)) {
        result <- list(success = FALSE, agent = key, message = worker$launch_error, rows = 0L)
      } else {
        alive <- tryCatch(worker$process$is_alive(), error = function(error) FALSE)
        if (alive && !file.exists(worker$result_file)) next
        result <- read_real_worker_result(worker)
        if (is.null(result)) next
      }

      worker$handled <- TRUE
      worker$result <- result
      workers[[key]] <- worker
      results[[key]] <- result
      changed <- TRUE
      if (isTRUE(result$success)) {
        send_progress(key, "completed", paste(format(result$rows, big.mark = ","), "results"))
      } else {
        send_progress(key, "error", result$message)
      }
    }

    if (changed) {
      workers_value(workers)
      results_value(results)
    }
    if (length(workers) && all(vapply(workers, function(x) isTRUE(x$handled), logical(1)))) {
      running(FALSE)
      run_finished(Sys.time())
      showNotification("Analysis run finished.", type = "message")
    }
  })

  output$run_overview <- renderUI({
    results <- results_value()
    total <- length(selected_value())
    completed <- sum(vapply(results, function(x) isTRUE(x$success), logical(1)))
    failed <- length(results) - completed
    elapsed <- if (is.null(run_started())) 0 else as.integer(difftime(run_finished() %||% Sys.time(), run_started(), units = "secs"))
    div(
      class = "run-overview",
      span(class = "eyebrow", if (running()) "PIPELINE RUNNING" else if (total) "RUN COMPLETE" else "READY"),
      div(class = "overview-number", paste0(completed, "/", total)),
      p("modules completed successfully"),
      div(class = "overview-stats", span(paste(elapsed, "sec")), span(paste(failed, "failed")))
    )
  })

  `%||%` <- function(value, fallback) if (is.null(value) || !length(value)) fallback else value

  result_table <- function(result) {
    if (!isTRUE(result$success)) return(NULL)
    table <- result$result$results
    if (is.matrix(table)) table <- data.frame(Pathway = rownames(table), table, check.names = FALSE)
    if (!is.data.frame(table)) return(NULL)
    table
  }

  output$results_section <- renderUI({
    results <- results_value()
    if (!length(results)) return(NULL)
    cards <- lapply(names(results), function(key) {
      result <- results[[key]]
      data <- result_table(result)
      meta <- module_meta[[key]]
      div(
        class = paste("result-card", if (!isTRUE(result$success)) "result-card-error"),
        div(class = "result-card-head", div(span(class = "result-kicker", meta$input), h3(meta$title)), span(class = "result-count", if (isTRUE(result$success)) paste(result$rows, "results") else "Failed")),
        p(class = "result-message", result$message),
        if (isTRUE(result$success) && !is.null(data) && nrow(data)) {
          DTOutput(paste0("table_", key))
        } else if (isTRUE(result$success)) {
          div(class = "empty-result", "Analysis completed with no significant results.")
        } else {
          div(class = "empty-result", "The module did not complete. Review the input requirements and message above.")
        }
      )
    })
    section(class = "results-section", div(class = "section-heading", span(class = "step", "04"), div(h2("Research results"), p("Top findings from the completed scientific modules"))), div(class = "results-grid", cards))
  })

  observe({
    results <- results_value()
    for (key in names(results)) local({
      local_key <- key
      output[[paste0("table_", local_key)]] <- renderDT({
        data <- result_table(results_value()[[local_key]])
        validate(need(!is.null(data) && nrow(data), "No result rows"))
        DT::datatable(utils::head(data, 10), rownames = FALSE, options = list(dom = "t", scrollX = TRUE, pageLength = 10), class = "compact")
      })
    })
  })

  build_report <- function(file) {
    build_combined_html_report(
      destination = file,
      interpretation_bundle = results_controller$interpretation_bundle(),
      selected_agents = selected_value()
    )
  }

  output$download_report <- downloadHandler(
    filename = function() paste0("OncoProfilingTools_Report_", format(Sys.Date(), "%Y%m%d"), ".html"),
    content = build_report,
    contentType = "text/html"
  )

  output$final_summary <- renderUI({
    if (is.null(run_finished())) return(NULL)
    results <- results_value()
    selected <- selected_value()
    status_items <- lapply(names(module_meta), function(key) {
      state <- if (key %in% names(results) && isTRUE(results[[key]]$success)) "Completed" else if (key %in% selected) "Failed" else "Not executed"
      div(class = paste("summary-module", tolower(gsub(" ", "-", state))), span(if (state == "Completed") "✓" else if (state == "Failed") "!" else "–"), module_meta[[key]]$title, tags$small(state))
    })
    section(
      class = "panel final-summary",
      div(class = "summary-head", div(span(class = "eyebrow", "FINAL SUMMARY"), h2("Analysis run complete"), p(input$dataset$name)), downloadButton("download_report", "Download report", class = "report-button")),
      div(class = "summary-metrics", div(span("Dataset"), strong(input$dataset$name)), div(span("Analysis genes"), strong(paste0(format(length(profile()$analysis_genes), big.mark = ","), " / ", format(length(profile()$genes), big.mark = ",")))), div(span("Modules executed"), strong(length(selected)))),
      div(class = "summary-modules", status_items)
    )
  })
}

shinyApp(ui, server)
