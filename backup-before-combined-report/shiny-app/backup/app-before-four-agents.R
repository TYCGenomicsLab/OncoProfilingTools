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
      href = "styles.css"
    )
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
