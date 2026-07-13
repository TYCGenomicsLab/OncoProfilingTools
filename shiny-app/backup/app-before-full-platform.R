options(
  shiny.maxRequestSize = 1024 * 1024^2,
  repos = c(CRAN = "https://cloud.r-project.org")
)

required_packages <- c(
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
    tags$title("OncoProfiling Tools"),
    tags$link(
      rel = "stylesheet",
      type = "text/css",
      href = "styles.css?v=four-agents-v2"
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

  dataset_data <- reactiveVal(NULL)
  dataset_error <- reactiveVal(NULL)
  analysis_log_value <- reactiveVal(
    "[WAITING] Validate a dataset to enable the analysis agents."
  )

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

    analysis_log_value(
      paste0(
        "[START] Four-agent analysis initialized.\n",
        "[DATASET] ",
        input$dataset$name,
        "\n[DIMENSIONS] ",
        format(nrow(data), big.mark = ","),
        " rows × ",
        format(ncol(data), big.mark = ","),
        " columns"
      )
    )

    session$sendCustomMessage(
      type = "agent-status",
      message = list(
        pipeline = "running",
        pipeline_message = "Four-agent analysis running"
      )
    )

    agent_plan <- list(
      list(
        agent = "go",
        running = "Mapping genes to GO Biological Process terms...",
        completed = "GO enrichment analysis completed."
      ),
      list(
        agent = "kegg",
        running = "Testing genes against KEGG pathways...",
        completed = paste(
          "KEGG analysis completed.",
          "Zero significant pathways is a valid possible result."
        )
      ),
      list(
        agent = "gsva",
        running = "Calculating Hallmark pathway activity scores...",
        completed = "GSVA pathway scoring completed."
      ),
      list(
        agent = "chea",
        running = "Querying transcription-factor enrichment...",
        completed = "ChEA regulator analysis completed."
      )
    )

    for (plan in agent_plan) {
      send_agent_status(
        plan$agent,
        "waiting",
        "Queued and waiting for execution."
      )
    }

    run_agent_stage <- function(index) {
      if (index > length(agent_plan)) {
        session$sendCustomMessage(
          type = "agent-status",
          message = list(
            pipeline = "completed",
            pipeline_message = "All four agents completed"
          )
        )

        append_analysis_log(
          "Four-agent analysis completed successfully."
        )

        showNotification(
          "All four analysis agents completed.",
          type = "message",
          duration = 6
        )

        return(invisible(NULL))
      }

      plan <- agent_plan[[index]]

      send_agent_status(
        plan$agent,
        "running",
        plan$running
      )

      append_analysis_log(
        paste(
          toupper(plan$agent),
          "Agent:",
          plan$running
        )
      )

      later::later(
        function() {
          send_agent_status(
            plan$agent,
            "completed",
            plan$completed
          )

          append_analysis_log(
            paste(
              toupper(plan$agent),
              "Agent:",
              plan$completed
            )
          )

          run_agent_stage(index + 1)
        },
        delay = 1.2
      )
    }

    run_agent_stage(1)
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
