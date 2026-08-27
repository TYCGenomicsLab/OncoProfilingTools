options(shiny.maxRequestSize = 1024 * 1024^2)

suppressPackageStartupMessages({
  library(shiny)
  library(readr)
  library(DT)
})

source("real_pipeline.R", local = TRUE)
source("run_real_agents.R", local = TRUE)
source("interpretation.R", local = TRUE)
source("interactive_visuals.R", local = TRUE)
source("results_helpers.R", local = TRUE)
source("production_reporting.R", local = TRUE)

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
  read_analysis_dataset(path, name)
}

dataset_profile <- function(data, pvalue_cutoff = 0.05, effect_cutoff = 1, max_genes = 2000L) {
  empty <- list(
    gene = FALSE, expression = FALSE, drug = FALSE,
    genes = character(), analysis_genes = character(), gene_column = NULL,
    selection_pvalue_column = NULL, selection_effect_column = NULL,
    original_gene_count = NULL,
    selection_note = NULL, selection_error = NULL, mapping_metadata = NULL,
    mapping = data.frame(), import_metadata = attr(data, "import_metadata")
  )
  if (is.null(data) || nrow(data) == 0L || ncol(data) == 0L) return(empty)

  normalized <- normalise_column_name(names(data))
  gene_column <- detect_gene_column(data)
  genes <- if (!is.null(gene_column)) unique(trimws(as.character(data[[gene_column]]))) else character()
  genes <- genes[!is.na(genes) & nzchar(genes)]
  recognized_gene_headers <- c(
    "genesymbol", "hugosymbol", "gene", "genes", "symbol", "genename",
    "geneid", "ensembl", "ensemblgeneid", "entrezid", "hgncsymbol",
    "externalgenename"
  )
  explicit_gene_column <- !is.null(gene_column) &&
    normalise_column_name(gene_column) %in% recognized_gene_headers

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
  if (length(genes) >= 2L && isTRUE(explicit_gene_column)) {
    selection <- tryCatch(
      prepare_gene_input(
        data,
        pvalue_cutoff = pvalue_cutoff,
        effect_cutoff = effect_cutoff,
        max_genes = max_genes
      ),
      error = function(error) {
        selection_error <<- conditionMessage(error)
        NULL
      }
    )
  }
  analysis_genes <- if (is.null(selection)) character() else selection$genes

  list(
    gene = length(analysis_genes) >= 2L && isTRUE(explicit_gene_column),
    expression = expression,
    drug = (compound && response) || wide_prism,
    genes = genes,
    analysis_genes = analysis_genes,
    gene_column = gene_column,
    selection_pvalue_column = if (is.null(selection)) NULL else selection$pvalue_column,
    selection_effect_column = if (is.null(selection)) NULL else selection$effect_column,
    original_gene_count = if (is.null(selection)) NULL else selection$original_gene_count,
    selection_note = if (is.null(selection)) NULL else selection$selection_note,
    selection_error = selection_error,
    mapping_metadata = if (is.null(selection)) NULL else selection$mapping_metadata,
    mapping = if (is.null(selection)) data.frame() else selection$mapping,
    import_metadata = attr(data, "import_metadata")
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

source("workflow_ui.R", local = TRUE)

ui <- build_app_ui()

server <- function(input, output, session) {
  data_value <- reactiveVal(NULL)
  data_error <- reactiveVal(NULL)
  workers_value <- reactiveVal(list())
  results_value <- reactiveVal(list())
  selected_value <- reactiveVal(character())
  run_started <- reactiveVal(NULL)
  run_finished <- reactiveVal(NULL)
  running <- reactiveVal(FALSE)
  workflow_mode <- reactiveVal("home")
  progress_states <- reactiveVal(setNames(vector("list", length(module_meta)), names(module_meta)))

  route_agents <- reactive({
    if (identical(workflow_mode(), "drug")) drug_agent_keys else biomarker_agent_keys
  })

  ollama_settings <- reactive({
    provider <- input$interpretation_provider %or_else% "ollama"
    list(
      enabled = !identical(provider, "computed"),
      provider = provider,
      host = input$ollama_host %or_else% "http://127.0.0.1:11434",
      model = input$ollama_model %or_else% "llama3.1:8b",
      timeout_seconds = input$ollama_timeout %or_else% 300,
      num_predict = 4096,
      openai_model = input$openai_model %or_else% "gpt-5.6-terra",
      openai_reasoning_effort = input$openai_reasoning %or_else% "medium",
      openai_timeout_seconds = input$openai_timeout %or_else% 240,
      openai_max_output_tokens = 12000L,
      openai_data_consent = isTRUE(input$openai_data_consent)
    )
  })

  output$openai_key_status <- renderUI({
    configured <- nzchar(trimws(Sys.getenv("OPENAI_API_KEY", "")))
    span(
      class = paste("openai-key-chip", if (configured) "openai-key-ready" else "openai-key-missing"),
      if (configured) "API key detected" else "API key not configured"
    )
  })

  output$ai_privacy_status <- renderUI({
    provider <- input$interpretation_provider %or_else% "ollama"
    if (provider %in% c("openai", "compare")) {
      div(class = "system-status system-status-external", span(), if (isTRUE(provider == "compare")) "Local + external comparison" else "External AI enabled")
    } else {
      div(class = "system-status", span(), "Local-only workspace")
    }
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

  output$app_content <- renderUI({
    if (identical(workflow_mode(), "home")) landing_page_ui() else workflow_page_ui(workflow_mode())
  })

  navigate_to <- function(mode) {
    reset_analysis_state(clear_files = TRUE)
    data_value(NULL)
    data_error(NULL)
    workflow_mode(mode)
    updateQueryString(
      if (identical(mode, "home")) "?" else paste0("?workflow=", mode),
      mode = "replace",
      session = session
    )
  }

  observeEvent(input$nav_biomarker, navigate_to("biomarker"), ignoreInit = TRUE)
  observeEvent(input$nav_drug, navigate_to("drug"), ignoreInit = TRUE)
  observeEvent(input$back_home, navigate_to("home"), ignoreInit = TRUE)
  observeEvent(input$results_back_home, navigate_to("home"), ignoreInit = TRUE)
  observeEvent(session$clientData$url_search, {
    query <- session$clientData$url_search %or_else% ""
    if (grepl("workflow=drug", query, fixed = TRUE)) workflow_mode("drug")
    if (grepl("workflow=biomarker", query, fixed = TRUE)) workflow_mode("biomarker")
  }, once = TRUE)

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
    progress_states(setNames(vector("list", length(module_meta)), names(module_meta)))
    if (isTRUE(clear_files)) clear_previous_artifacts()
    session$sendCustomMessage("reset-analysis-ui", list())
    invisible(NULL)
  }

  # Never expose output files left behind by an earlier app session.
  clear_previous_artifacts()

  session$onSessionEnded(function() {
    stop_active_workers()
  })

  output$biomarker_configuration <- renderUI({
    data <- data_value()
    if (is.null(data)) {
      return(div(class = "scientific-note", "Upload data to reveal relevant DEG filters."))
    }
    normalized <- normalise_column_name(names(data))
    direction_index <- which(normalized %in% c("nominalsign", "direction", "regulation", "group"))
    type_index <- which(normalized %in% c("genetype", "biotype"))
    direction_control <- NULL
    if (length(direction_index)) {
      values <- sort(unique(as.character(data[[direction_index[[1L]]]])))
      values <- values[!is.na(values) & nzchar(values)]
      default <- if ("higherSMI" %in% values) "higherSMI" else "__all__"
      direction_control <- selectInput(
        "gene_direction", "Gene group",
        choices = c("All values" = "__all__", values), selected = default
      )
    }
    type_control <- NULL
    if (length(type_index)) {
      types <- unique(tolower(as.character(data[[type_index[[1L]]]])))
      type_control <- checkboxInput(
        "protein_coding_only", "Protein-coding genes only",
        value = "protein_coding" %in% types
      )
    }
    div(
      class = "focused-config",
      div(class = "focused-config-title", strong("Focused DEG configuration"), span("Applied before enrichment")),
      div(
        class = "focused-config-grid",
        direction_control,
        type_control,
        numericInput("pvalue_cutoff", "Adjusted p/FDR cutoff", value = 0.05, min = 0.0001, max = 1, step = 0.01),
        numericInput("effect_cutoff", "Absolute effect cutoff", value = 1, min = 0, max = 20, step = 0.25),
        numericInput("max_genes", "Maximum genes", value = 2000, min = 2, max = 10000, step = 100),
        textAreaInput(
          "experimental_design",
          "Experimental design / biological comparison (optional)",
          value = "",
          placeholder = "Example: cases versus matched controls; tissue, time point, and covariates",
          rows = 2
        )
      )
    )
  })

  configured_data <- reactive({
    data <- data_value()
    if (is.null(data) || !identical(workflow_mode(), "biomarker")) return(data)
    normalized <- normalise_column_name(names(data))
    keep <- rep(TRUE, nrow(data))
    direction_index <- which(normalized %in% c("nominalsign", "direction", "regulation", "group"))
    direction <- input$gene_direction %or_else% "__all__"
    if (length(direction_index) && !identical(direction, "__all__")) {
      keep <- keep & as.character(data[[direction_index[[1L]]]]) == direction
    }
    type_index <- which(normalized %in% c("genetype", "biotype"))
    if (length(type_index) && isTRUE(input$protein_coding_only)) {
      keep <- keep & tolower(as.character(data[[type_index[[1L]]]])) == "protein_coding"
    }
    filtered <- data[keep, , drop = FALSE]
    attr(filtered, "import_metadata") <- attr(data, "import_metadata")
    filtered
  })

  profile <- reactive(dataset_profile(
    configured_data(),
    pvalue_cutoff = input$pvalue_cutoff %or_else% 0.05,
    effect_cutoff = input$effect_cutoff %or_else% 1,
    max_genes = input$max_genes %or_else% 2000L
  ))

  compatibility <- reactive({
    p <- profile()
    compatible <- c(
      go = p$gene, kegg = p$gene, reactome = p$gene, wikipathways = p$gene, string = p$gene, hallmark = p$gene, chea = p$gene,
      gsva = p$expression, immune = p$expression, drug = p$drug
    )
    compatible[names(compatible) %in% route_agents()]
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
    for (key in route_agents()) {
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
    import <- profile()$import_metadata
    header_note <- if (!is.null(import) && !isTRUE(import$header_detected)) " · headerless list detected" else ""
    div(class = "dataset-alert dataset-ready", span("✓"), strong(input$dataset$name), paste0(" validated", header_note))
  })

  output$dataset_metrics <- renderUI({
    data <- configured_data()
    if (is.null(data)) return(NULL)
    p <- profile()
    div(
      class = "dataset-metrics",
      div(span("Rows used"), strong(format(nrow(data), big.mark = ","))),
      div(span("Columns"), strong(format(ncol(data), big.mark = ","))),
      if (identical(workflow_mode(), "biomarker") && !is.null(p$mapping_metadata)) div(
        span("Mapped symbols"),
        strong(format(length(p$analysis_genes), big.mark = ","))
      ) else if (identical(workflow_mode(), "biomarker") && isTRUE(p$expression)) div(
        span("Matrix orientation"), strong("Samples × genes")
      ) else NULL,
      div(span("Compatible agents"), strong(sum(compatibility())))
    )
  })

  output$gene_selection_note <- renderUI({
    if (is.null(data_value()) || !identical(workflow_mode(), "biomarker")) return(NULL)
    p <- profile()
    if (!is.null(p$selection_error)) {
      return(div(class = "gene-selection-note gene-selection-warning", strong("Gene-list input needs attention"), p$selection_error))
    }
    if (is.null(p$selection_note)) return(NULL)
    div(
      class = "gene-selection-note",
      strong("Input and mapping summary"),
      p$selection_note,
      tags$small("For directional biology, upload up- and down-regulated signatures as separate runs.")
    )
  })

  send_progress <- function(key, status, message) {
    states <- progress_states()
    states[[key]] <- list(status = status, message = message)
    progress_states(states)
    session$sendCustomMessage("analysis-progress", list(key = key, status = status, message = message))
  }

  observeEvent(input$run_analysis, {
    data <- configured_data()
    if (is.null(data)) {
      showNotification("Upload a valid dataset first.", type = "warning")
      return()
    }
    if (running()) return()

    requested <- vapply(route_agents(), function(key) isTRUE(input[[paste0("module_", key)]]), logical(1))
    selected <- names(requested)[requested & compatibility()[names(requested)]]
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

    for (key in selected) send_progress(key, "queued", "Queued")

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
        start_real_agent_worker(
          key,
          analysis_run,
          pvalue_cutoff = input$pvalue_cutoff %or_else% 0.05,
          effect_cutoff = input$effect_cutoff %or_else% 1,
          max_genes = input$max_genes %or_else% 2000L
        ),
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

  output$run_status_bar <- renderUI({
    if (is.null(run_started())) return(NULL)
    selected <- selected_value()
    states <- progress_states()
    results <- results_value()
    completed <- sum(vapply(results, function(result) isTRUE(result$success), logical(1)))
    failed <- length(results) - completed
    elapsed <- as.integer(difftime(run_finished() %or_else% Sys.time(), run_started(), units = "secs"))
    rows <- lapply(selected, function(key) {
      state <- states[[key]] %or_else% list(status = "queued", message = "Queued")
      div(
        id = paste0("progress-", key),
        class = paste("progress-row", paste0("progress-", state$status)),
        span(
          class = "progress-indicator",
          if (state$status == "completed") "✓" else if (state$status == "error") "!" else if (state$status == "running") "↻" else "·"
        ),
        div(strong(module_meta[[key]]$title), span(class = "progress-message", state$message))
      )
    })
    section(
      class = "panel compact-run-status",
      div(
        class = "compact-run-summary",
        div(span(class = "eyebrow", if (running()) "ANALYSIS RUNNING" else "ANALYSIS COMPLETE"), strong(paste0(completed, "/", length(selected), " agents"))),
        div(span(paste(elapsed, "sec")), span(paste(failed, "failed")))
      ),
      div(class = "compact-progress-list", rows)
    )
  })

  run_context <- reactive({
    p <- profile()
    data <- configured_data()
    list(
      workflow = workflow_mode(),
      input = list(
        name = if (is.null(input$dataset)) NULL else input$dataset$name,
        size_bytes = if (is.null(input$dataset)) NULL else input$dataset$size,
        checksum_md5 = if (is.null(input$dataset) || !file.exists(input$dataset$datapath)) NULL else unname(tools::md5sum(input$dataset$datapath)),
        rows = if (is.null(data)) 0L else nrow(data),
        columns = if (is.null(data)) 0L else ncol(data),
        import = p$import_metadata
      ),
      mapping = p$mapping_metadata,
      mapping_table = p$mapping,
      configuration = list(
        gene_direction = input$gene_direction %or_else% "all",
        protein_coding_only = isTRUE(input$protein_coding_only),
        pvalue_cutoff = input$pvalue_cutoff %or_else% 0.05,
        effect_cutoff = input$effect_cutoff %or_else% 1,
        max_genes = input$max_genes %or_else% 2000L,
        gene_column = p$gene_column,
        pvalue_column = p$selection_pvalue_column,
        effect_column = p$selection_effect_column,
        original_gene_count = p$original_gene_count,
        selection_note = p$selection_note,
        experimental_design = input$experimental_design %or_else% "Not supplied",
        tested_gene_universe = "Not supplied; agent/package defaults were used."
      ),
      started_at = run_started(),
      finished_at = run_finished()
    )
  })

  output$results_center <- renderUI({
    if (is.null(run_started())) return(NULL)
    results_center_ui(selected_value(), workflow = workflow_mode())
  })

  results_controller <- register_results_server(
    input,
    output,
    session,
    active_agents = selected_value,
    analysis_complete = analysis_complete,
    ollama_settings = ollama_settings,
    run_context = run_context
  )

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
      selected_agents = selected_value(),
      run_context = run_context()
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
