options(shiny.maxRequestSize = 500 * 1024^2)

options(
  shiny.maxRequestSize = 1024 * 1024^2,
  repos = c(CRAN = "https://cloud.r-project.org")
)

required_packages <- c(
  "base64enc",
  "shiny",
  "readr",
  "dplyr",
  "DT"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(shiny)
library(readr)
library(dplyr)
library(DT)
source("results_helpers.R", local = TRUE)
source("real_pipeline.R", local = TRUE)
source("agent_inputs.R", local = TRUE)

read_uploaded_dataset <- function(file_path, file_name) {
  extension <- tolower(tools::file_ext(file_name))

  if (extension == "csv") {
    return(
      readr::read_csv(
        file_path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "unique"
      )
    )
  }

  if (extension %in% c("tsv", "txt")) {
    return(
      readr::read_tsv(
        file_path,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "unique"
      )
    )
  }

  stop("Unsupported file type. Upload CSV, TSV, or TXT.")
}

ui <- fluidPage(
  tags$head(
    tags$script(src = "custom-upload-fix.js"),
    tags$title("OncoProfiling Tools"),
    tags$link(
      rel = "stylesheet",
      type = "text/css",
      href = "styles.css?v=agent-input-controls-v1"
    ),
    tags$script(src = "status.js")
  ),

  div(
    class = "app-page",

    div(
      class = "top-header",

      div(
        class = "brand-area",
        div(class = "brand-logo", "OP"),
        div(
          h1("OncoProfiling Tools"),
          p("Multi-agent genomic interpretation workspace")
        )
      ),

      div(
        class = "header-badge",
        span(class = "header-dot"),
        "System Ready"
      )
    ),

    div(
      class = "live-metrics-grid",

      div(
        class = "live-metric-card",
        div(class = "live-metric-icon metric-dataset-icon", "✓"),
        div(
          class = "live-metric-content",
          span(class = "live-metric-label", "Dataset"),
          uiOutput("metric_dataset_status")
        )
      ),

      div(
        class = "live-metric-card",
        div(class = "live-metric-icon metric-gene-icon", "DNA"),
        div(
          class = "live-metric-content",
          span(class = "live-metric-label", "Genes / Features"),
          uiOutput("metric_gene_count")
        )
      ),

      div(
        class = "live-metric-card",
        div(class = "live-metric-icon metric-sample-icon", "S"),
        div(
          class = "live-metric-content",
          span(class = "live-metric-label", "Samples"),
          uiOutput("metric_sample_count")
        )
      ),

      div(
        class = "live-metric-card",
        div(class = "live-metric-icon metric-agent-icon", "AI"),
        div(
          class = "live-metric-content",
          span(class = "live-metric-label", "Agents"),
          uiOutput("metric_agent_count")
        )
      ),

      div(
        class = "live-metric-card",
        div(class = "live-metric-icon metric-runtime-icon", "T"),
        div(
          class = "live-metric-content",
          span(class = "live-metric-label", "Runtime"),
          uiOutput("metric_runtime")
        )
      )
    ),

    div(
      class = "main-grid",

      div(
        class = "upload-card glass-card",

        div(
          class = "section-label",
          "STEP 1"
        ),

        h2("Upload genomic dataset"),

        p(
          class = "section-description",
          paste(
            "Upload a CSV, TSV, or TXT dataset.",
            "The application will validate and preview it before analysis."
          )
        ),

        div(
          class = "upload-box",

          div(
            class = "upload-copy",
            div(class = "upload-icon", "⇧"),
            h3("Choose your dataset"),
            p("Supported formats: CSV, TSV, TXT")
          ),

          fileInput(
            inputId = "dataset",
            label = NULL,
            accept = c(".csv", ".tsv", ".txt")
          )
        ),

        uiOutput("uploaded_file_details"),

        actionButton(
          inputId = "validate_dataset",
          label = "Validate Dataset",
          class = "run-button"
        )
      ),

      div(
        class = "status-card glass-card",

        div(
          class = "section-label",
          "DATASET STATUS"
        ),

        h2("Validation summary"),

        uiOutput("validation_summary"),

        div(
          class = "status-log",
          verbatimTextOutput("status")
        )
      )
    ),

    agent_input_control_ui(),

    div(
      class = "agents-workspace glass-card",

      div(
        class = "agents-heading",

        div(
          div(class = "section-label", "STEP 2"),
          h2("Four intelligent analysis agents"),
          p(
            class = "section-description",
            paste(
              "Each agent runs a specialized genomic analysis.",
              "Live status, results, and interpretation are shown below."
            )
          )
        ),

        div(
          class = "pipeline-status",
          span(class = "pipeline-status-dot"),
          span(id = "pipeline-status-text", "Waiting for dataset")
        )
      ),

      div(
        class = "agents-grid",

        div(
          id = "go-agent-card",
          class = "agent-card agent-waiting",

          div(
            class = "agent-card-top",
            div(class = "agent-icon agent-icon-go", "GO"),
            span(
              id = "go-agent-badge",
              class = "agent-status-badge",
              "Waiting"
            )
          ),

          h3("GO Agent"),
          p("Identifies enriched Gene Ontology biological processes."),

          div(
            class = "agent-progress-track",
            div(
              id = "go-agent-progress",
              class = "agent-progress-fill"
            )
          ),

          div(
            id = "go-agent-message",
            class = "agent-message",
            "Waiting for analysis to begin."
          )
        ),

        div(
          id = "kegg-agent-card",
          class = "agent-card agent-waiting",

          div(
            class = "agent-card-top",
            div(class = "agent-icon agent-icon-kegg", "KG"),
            span(
              id = "kegg-agent-badge",
              class = "agent-status-badge",
              "Waiting"
            )
          ),

          h3("KEGG Agent"),
          p("Detects significantly enriched molecular pathways."),

          div(
            class = "agent-progress-track",
            div(
              id = "kegg-agent-progress",
              class = "agent-progress-fill"
            )
          ),

          div(
            id = "kegg-agent-message",
            class = "agent-message",
            "Waiting for analysis to begin."
          )
        ),

        div(
          id = "gsva-agent-card",
          class = "agent-card agent-waiting",

          div(
            class = "agent-card-top",
            div(class = "agent-icon agent-icon-gsva", "GS"),
            span(
              id = "gsva-agent-badge",
              class = "agent-status-badge",
              "Waiting"
            )
          ),

          h3("GSVA Agent"),
          p("Calculates sample-level Hallmark pathway activity scores."),

          div(
            class = "agent-progress-track",
            div(
              id = "gsva-agent-progress",
              class = "agent-progress-fill"
            )
          ),

          div(
            id = "gsva-agent-message",
            class = "agent-message",
            "Waiting for analysis to begin."
          )
        ),

        div(
          id = "chea-agent-card",
          class = "agent-card agent-waiting",

          div(
            class = "agent-card-top",
            div(class = "agent-icon agent-icon-chea", "TF"),
            span(
              id = "chea-agent-badge",
              class = "agent-status-badge",
              "Waiting"
            )
          ),

          h3("ChEA Agent"),
          p("Discovers candidate transcription-factor regulators."),

          div(
            class = "agent-progress-track",
            div(
              id = "chea-agent-progress",
              class = "agent-progress-fill"
            )
          ),

          div(
            id = "chea-agent-message",
            class = "agent-message",
            "Waiting for analysis to begin."
          )
        )
      ),

      actionButton(
        inputId = "run_analysis",
        label = "Run Four-Agent Analysis",
        class = "run-button run-all-agents-button"
      ),

      div(
        class = "analysis-log-panel",
        div(
          class = "analysis-log-heading",
          span("LIVE PIPELINE LOG"),
          span(class = "live-log-indicator")
        ),
        verbatimTextOutput("analysis_log")
      )
    ),

    results_center_ui(),

    div(
      class = "preview-card glass-card",

      div(
        class = "preview-heading",

        div(
          div(class = "section-label", "DATA PREVIEW"),
          h2("Uploaded dataset")
        ),

        uiOutput("dataset_dimensions")
      ),

      DTOutput("dataset_preview")
    )
  )
)

server <- function(input, output, session) {

  register_results_server(output, session)


  dataset_data <- reactiveVal(NULL)
  dataset_error <- reactiveVal(NULL)

  dataset_profile <- reactive({
    detect_agent_compatibility(dataset_data())
  })
  completed_agents <- reactiveVal(0L)
  analysis_start_time <- reactiveVal(NULL)
  analysis_runtime_seconds <- reactiveVal(0L)

  analysis_log_value <- reactiveVal(
    "[WAITING] Validate a dataset to enable the analysis agents."
  )
  real_workers <- reactiveVal(list())
  real_pipeline_running <- reactiveVal(FALSE)
  real_pipeline_finished <- reactiveVal(FALSE)


  send_agent_status <- function(agent, status, message) {
    session$sendCustomMessage(
      type = "agent-status",
      message = list(
        agent = agent,
        status = status,
        message = message
      )
    )
  }

  append_analysis_log <- function(message) {
    current_log <- analysis_log_value()

    timestamp <- format(
      Sys.time(),
      "%H:%M:%S"
    )

    analysis_log_value(
      paste0(
        current_log,
        "\n[",
        timestamp,
        "] ",
        message
      )
    )
  }

  output$analysis_log <- renderText({
    analysis_log_value()
  })

  output$metric_dataset_status <- renderUI({
    if (is.null(dataset_data())) {
      return(
        div(
          class = "live-metric-value metric-value-waiting",
          "Waiting"
        )
      )
    }

    div(
      class = "live-metric-value metric-value-ready",
      "Validated"
    )
  })

  output$metric_gene_count <- renderUI({
    data <- dataset_data()

    value <- if (is.null(data)) {
      "—"
    } else {
      format(ncol(data), big.mark = ",")
    }

    div(class = "live-metric-value", value)
  })

  output$metric_sample_count <- renderUI({
    data <- dataset_data()

    value <- if (is.null(data)) {
      "—"
    } else {
      format(nrow(data), big.mark = ",")
    }

    div(class = "live-metric-value", value)
  })

  output$metric_agent_count <- renderUI({
    div(
      class = "live-metric-value",
      paste0(completed_agents(), " / 4")
    )
  })

  output$metric_runtime <- renderUI({
    seconds <- analysis_runtime_seconds()

    minutes <- seconds %/% 60
    remaining_seconds <- seconds %% 60

    div(
      class = "live-metric-value",
      sprintf("%02d:%02d", minutes, remaining_seconds)
    )
  })

  output$detected_dataset_summary <- renderUI({
    profile <- dataset_profile()

    div(
      class = "detected-dataset-summary",

      div(
        class = "detected-summary-item",
        span("Detected type"),
        strong(profile$dataset_type)
      ),

      div(
        class = "detected-summary-item",
        span("Gene identifier"),
        strong(
          if (is.null(profile$gene_column)) {
            "Not detected"
          } else {
            profile$gene_column
          }
        )
      ),

      div(
        class = "detected-summary-item",
        span("Numeric samples"),
        strong(profile$sample_count)
      )
    )
  })

  output$selected_agent_summary <- renderUI({
    selected <- c(
      go = isTRUE(input$enable_go),
      kegg = isTRUE(input$enable_kegg),
      gsva = isTRUE(input$enable_gsva),
      chea = isTRUE(input$enable_chea)
    )

    compatible <- dataset_profile()$agent_compatibility

    valid_selected <- selected & compatible

    div(
      class = "selected-agent-summary",
      span(class = "selected-agent-dot"),
      strong(sum(valid_selected)),
      span("compatible agents selected")
    )
  })

  for (agent_key in names(agent_requirements)) {
    local({
      key <- agent_key

      output[[paste0(key, "_compatibility")]] <- renderUI({
        profile <- dataset_profile()
        compatible <- isTRUE(
          profile$agent_compatibility[[key]]
        )

        div(
          class = paste(
            "agent-compatibility",
            if (compatible) {
              "agent-compatible"
            } else {
              "agent-incompatible"
            }
          ),

          span(
            class = "agent-compatibility-dot"
          ),

          div(
            strong(
              if (compatible) {
                "Compatible"
              } else {
                "Not compatible"
              }
            ),
            span(profile$messages[[key]])
          )
        )
      })
    })
  }

  observe({
    profile <- dataset_profile()

    for (agent_key in names(agent_requirements)) {
      compatible <- isTRUE(
        profile$agent_compatibility[[agent_key]]
      )

      shinyjs_state <- compatible &&
        !is.null(dataset_data())

      session$sendCustomMessage(
        "agent-input-state",
        list(
          agent = agent_key,
          enabled = shinyjs_state
        )
      )
    }
  })

  output$uploaded_file_details <- renderUI({
    if (is.null(input$dataset)) {
      return(
        div(
          class = "file-placeholder",
          "No file selected yet."
        )
      )
    }

    div(
      class = "file-details",

      div(
        class = "file-type-icon",
        toupper(tools::file_ext(input$dataset$name))
      ),

      div(
        class = "file-copy",
        strong(input$dataset$name),
        span(
          paste(
            round(input$dataset$size / 1024, 2),
            "KB"
          )
        )
      )
    )
  })

  observeEvent(input$validate_dataset, {

    dataset_error(NULL)
    dataset_data(NULL)

    if (is.null(input$dataset)) {
      dataset_error("Upload a dataset before validation.")

      showNotification(
        "Please upload a dataset first.",
        type = "warning"
      )

      return()
    }

    tryCatch(
      {
        loaded_data <- read_uploaded_dataset(
          input$dataset$datapath,
          input$dataset$name
        )

        if (nrow(loaded_data) == 0) {
          stop("The uploaded dataset has no data rows.")
        }

        if (ncol(loaded_data) == 0) {
          stop("The uploaded dataset has no columns.")
        }

        dataset_data(loaded_data)

        showNotification(
          "Dataset validated successfully.",
          type = "message"
        )
      },
      error = function(error) {
        dataset_error(conditionMessage(error))

        showNotification(
          conditionMessage(error),
          type = "error",
          duration = 8
        )
      }
    )
  })

  observeEvent(input$run_analysis, {

    data <- dataset_data()

    if (is.null(data)) {
      showNotification(
        "Validate your dataset before running the agents.",
        type = "warning",
        duration = 6
      )

      return()
    }

    if (isTRUE(real_pipeline_running())) {
      showNotification(
        "An analysis is already running.",
        type = "warning",
        duration = 5
      )

      return()
    }

    profile <- dataset_profile()

    requested_agents <- c(
      go = isTRUE(input$enable_go),
      kegg = isTRUE(input$enable_kegg),
      gsva = isTRUE(input$enable_gsva),
      chea = isTRUE(input$enable_chea)
    )

    compatible_agents <- vapply(
      names(requested_agents),
      function(agent) {
        isTRUE(
          profile$agent_compatibility[[agent]]
        )
      },
      logical(1)
    )

    selected_agents <- names(
      requested_agents[
        requested_agents &
          compatible_agents
      ]
    )

    if (length(selected_agents) == 0) {
      showNotification(
        paste(
          "Select at least one compatible analysis agent",
          "before starting the pipeline."
        ),
        type = "warning",
        duration = 7
      )

      return()
    }

    completed_agents(0L)
    analysis_start_time(Sys.time())
    analysis_runtime_seconds(0L)

    real_pipeline_running(TRUE)
    real_pipeline_finished(FALSE)

    analysis_log_value(
      paste0(
        "[START] Real multi-agent analysis initialized.\n",
        "[DATASET] ",
        input$dataset$name,
        "\n[DIMENSIONS] ",
        format(nrow(data), big.mark = ","),
        " rows × ",
        format(ncol(data), big.mark = ","),
        " columns\n[AGENTS] ",
        paste(
          toupper(selected_agents),
          collapse = ", "
        )
      )
    )

    session$sendCustomMessage(
      type = "agent-status",
      message = list(
        pipeline = "running",
        pipeline_message = paste(
          length(selected_agents),
          "real analysis agents running"
        )
      )
    )

    all_agents <- c(
      "go",
      "kegg",
      "gsva",
      "chea"
    )

    for (agent in all_agents) {
      if (agent %in% selected_agents) {
        send_agent_status(
          agent,
          "waiting",
          "Preparing real analysis worker."
        )
      } else {
        send_agent_status(
          agent,
          "waiting",
          "Not selected for this analysis."
        )
      }
    }

    analysis_run <- tryCatch(
      create_real_analysis_run(data),
      error = function(error) {
        real_pipeline_running(FALSE)

        append_analysis_log(
          paste(
            "Unable to prepare the real analysis run:",
            conditionMessage(error)
          )
        )

        showNotification(
          conditionMessage(error),
          type = "error",
          duration = 10
        )

        return(NULL)
      }
    )

    if (is.null(analysis_run)) {
      return()
    }

    workers <- list()

    for (agent in selected_agents) {

      send_agent_status(
        agent,
        "running",
        paste(
          "Running real",
          toupper(agent),
          "analysis..."
        )
      )

      append_analysis_log(
        paste(
          toupper(agent),
          "worker started."
        )
      )

      worker <- tryCatch(
        start_real_agent_worker(
          agent = agent,
          analysis_run = analysis_run,
          pvalue_cutoff = 0.05
        ),
        error = function(error) {
          list(
            agent = agent,
            launch_error = conditionMessage(error),
            handled = FALSE
          )
        }
      )

      worker$handled <- FALSE
      workers[[agent]] <- worker
    }

    real_workers(workers)
  })


  observe({

    if (!isTRUE(real_pipeline_running())) {
      return()
    }

    invalidateLater(
      750,
      session
    )

    workers <- real_workers()

    if (length(workers) == 0) {
      return()
    }

    state_changed <- FALSE

    for (agent in names(workers)) {

      worker <- workers[[agent]]

      if (isTRUE(worker$handled)) {
        next
      }

      if (!is.null(worker$launch_error)) {

        send_agent_status(
          agent,
          "error",
          worker$launch_error
        )

        append_analysis_log(
          paste(
            toupper(agent),
            "worker could not start:",
            worker$launch_error
          )
        )

        worker$handled <- TRUE
        workers[[agent]] <- worker
        state_changed <- TRUE

        next
      }

      process_alive <- tryCatch(
        worker$process$is_alive(),
        error = function(error) {
          FALSE
        }
      )

      result_ready <- file.exists(
        worker$result_file
      )

      if (process_alive && !result_ready) {
        next
      }

      result <- read_real_worker_result(
        worker
      )

      if (is.null(result)) {
        next
      }

      if (isTRUE(result$success)) {

        completed_agents(
          completed_agents() + 1L
        )

        send_agent_status(
          agent,
          "completed",
          paste0(
            toupper(agent),
            " completed with ",
            format(
              result$rows,
              big.mark = ","
            ),
            " result rows."
          )
        )

        append_analysis_log(
          paste(
            toupper(agent),
            "Agent completed:",
            result$message
          )
        )

      } else {

        send_agent_status(
          agent,
          "error",
          result$message
        )

        append_analysis_log(
          paste(
            toupper(agent),
            "Agent failed:",
            result$message
          )
        )
      }

      worker$handled <- TRUE
      worker$result <- result
      workers[[agent]] <- worker
      state_changed <- TRUE
    }

    if (state_changed) {
      real_workers(workers)
    }

    all_finished <- length(workers) > 0 &&
      all(
        vapply(
          workers,
          function(worker) {
            isTRUE(worker$handled)
          },
          logical(1)
        )
      )

    if (!all_finished) {
      return()
    }

    real_pipeline_running(FALSE)
    real_pipeline_finished(TRUE)

    successful_agents <- sum(
      vapply(
        workers,
        function(worker) {
          !is.null(worker$result) &&
            isTRUE(worker$result$success)
        },
        logical(1)
      )
    )

    failed_agents <- length(workers) -
      successful_agents

    analysis_runtime_seconds(
      as.integer(
        difftime(
          Sys.time(),
          analysis_start_time(),
          units = "secs"
        )
      )
    )

    pipeline_message <- if (failed_agents == 0) {
      paste(
        "All",
        successful_agents,
        "selected agents completed"
      )
    } else {
      paste(
        successful_agents,
        "completed and",
        failed_agents,
        "failed"
      )
    }

    session$sendCustomMessage(
      type = "agent-status",
      message = list(
        pipeline = if (failed_agents == 0) {
          "completed"
        } else {
          "error"
        },
        pipeline_message = pipeline_message
      )
    )

    append_analysis_log(
      paste(
        "Real multi-agent analysis finished.",
        pipeline_message
      )
    )

    showNotification(
      pipeline_message,
      type = if (failed_agents == 0) {
        "message"
      } else {
        "warning"
      },
      duration = 8
    )
  })

  observe({
    invalidateLater(1000, session)

    started <- analysis_start_time()

    if (!is.null(started) && completed_agents() < 4L) {
      elapsed <- as.integer(
        difftime(
          Sys.time(),
          started,
          units = "secs"
        )
      )

      analysis_runtime_seconds(max(elapsed, 0L))
    }
  })

  output$status <- renderText({
    if (!is.null(dataset_error())) {
      return(
        paste(
          "[ERROR]",
          dataset_error()
        )
      )
    }

    data <- dataset_data()

    if (is.null(input$dataset)) {
      return("[WAITING] No dataset uploaded.")
    }

    if (is.null(data)) {
      return(
        paste(
          "[READY] File selected:",
          input$dataset$name,
          "\nClick Validate Dataset."
        )
      )
    }

    paste0(
      "[SUCCESS] Dataset loaded\n",
      "Rows: ", format(nrow(data), big.mark = ","), "\n",
      "Columns: ", format(ncol(data), big.mark = ","), "\n",
      "File: ", input$dataset$name
    )
  })

  output$validation_summary <- renderUI({
    data <- dataset_data()

    if (!is.null(dataset_error())) {
      return(
        div(
          class = "validation-state validation-error",
          div(class = "validation-state-icon", "×"),
          div(
            strong("Validation failed"),
            span(dataset_error())
          )
        )
      )
    }

    if (is.null(data)) {
      return(
        div(
          class = "validation-state validation-waiting",
          div(class = "validation-state-icon", "…"),
          div(
            strong("Waiting for validation"),
            span("Upload a file and click Validate Dataset.")
          )
        )
      )
    }

    div(
      class = "validation-state validation-success",
      div(class = "validation-state-icon", "✓"),
      div(
        strong("Dataset ready"),
        span(
          paste(
            format(nrow(data), big.mark = ","),
            "rows and",
            format(ncol(data), big.mark = ","),
            "columns detected."
          )
        )
      )
    )
  })

  output$dataset_dimensions <- renderUI({
    data <- dataset_data()

    if (is.null(data)) {
      return(
        span(
          class = "dimension-badge",
          "No data"
        )
      )
    }

    span(
      class = "dimension-badge dimension-ready",
      paste(
        format(nrow(data), big.mark = ","),
        "×",
        format(ncol(data), big.mark = ",")
      )
    )
  })


  # -----------------------------------------------------------------------
  # REAL AGENT VISUALIZATIONS
  # -----------------------------------------------------------------------

  project_root <- normalizePath(
    file.path(getwd(), ".."),
    mustWork = TRUE
  )

  agent_visualization_paths <- list(
    go = file.path(
      project_root,
      "output",
      "visualizations",
      "go_biological_process_dotplot.png"
    ),
    kegg = file.path(
      project_root,
      "output",
      "visualizations",
      "kegg_pathway_dotplot.png"
    ),
    gsva = file.path(
      project_root,
      "output",
      "gsva_bowel",
      "gsva_hallmark_heatmap.png"
    ),
    chea = file.path(
      project_root,
      "output",
      "chea_cms4",
      "chea_tf_dotplot.png"
    )
  )

  visualization_panel <- function(
      path,
      title,
      unavailable_message
  ) {
    if (!file.exists(path)) {
      return(
        div(
          class = "visualization-empty-state",
          div(class = "visualization-empty-icon", "⌁"),
          h4(title),
          p(unavailable_message)
        )
      )
    }

    extension <- tools::file_ext(path)

    temporary_plot <- tempfile(
      pattern = "agent-visualization-",
      fileext = paste0(".", extension)
    )

    file.copy(
      from = path,
      to = temporary_plot,
      overwrite = TRUE
    )

    image_data <- base64enc::dataURI(
      file = temporary_plot,
      mime = paste0("image/", extension)
    )

    tags$img(
      src = image_data,
      alt = title,
      class = "agent-result-image"
    )
  }

  output$go_visualization <- renderUI({
    visualization_panel(
      agent_visualization_paths$go,
      "GO visualization unavailable",
      paste(
        "Run the GO analysis and visualization pipeline",
        "to generate the biological-process dot plot."
      )
    )
  })

  output$kegg_visualization <- renderUI({
    visualization_panel(
      agent_visualization_paths$kegg,
      "No significant KEGG visualization",
      paste(
        "The current KEGG results contain zero enriched pathways.",
        "This is a valid biological result rather than an application error."
      )
    )
  })

  output$gsva_visualization <- renderUI({
    visualization_panel(
      agent_visualization_paths$gsva,
      "GSVA visualization unavailable",
      paste(
        "Run GSVA using a valid gene-expression matrix",
        "to generate the Hallmark pathway heatmap."
      )
    )
  })

  output$chea_visualization <- renderUI({
    visualization_panel(
      agent_visualization_paths$chea,
      "ChEA visualization unavailable",
      paste(
        "Run ChEA enrichment to generate",
        "the transcription-factor dot plot."
      )
    )
  })


  output$dataset_preview <- renderDT({
    data <- dataset_data()

    validate(
      need(
        !is.null(data),
        "Validate a dataset to see its preview."
      )
    )

    preview_data <- data |>
      head(20)

    datatable(
      preview_data,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "tip"
      ),
      class = "stripe hover compact"
    )
  })
}

shinyApp(ui = ui, server = server)
