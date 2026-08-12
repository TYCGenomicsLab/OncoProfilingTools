# OncoProfiling Tools REST API

suppressPackageStartupMessages({
  library(plumber)
  library(jsonlite)
})

api_directory <- normalizePath(
  ".",
  mustWork = TRUE
)

api_project_root <- normalizePath(
  file.path(api_directory, ".."),
  mustWork = TRUE
)

original_working_directory <- getwd()

setwd(api_project_root)

source(
  file.path(
    "shiny-app",
    "run_real_agents.R"
  ),
  local = TRUE
)

setwd(original_working_directory)

#* API information
#* @get /
function() {
  list(
    success = TRUE,
    application = "OncoProfiling Tools API",
    version = "1.0.0",
    agents = c(
      "go",
      "kegg",
      "gsva",
      "chea"
    )
  )
}

#* Backend health check
#* @get /health
function() {
  required_packages <- c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichR",
    "GSVA",
    "msigdbr",
    "pheatmap",
    "ggplot2",
    "readr"
  )

  package_status <- vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )

  list(
    success = all(package_status),
    status = if (all(package_status)) {
      "healthy"
    } else {
      "missing_packages"
    },
    packages = as.list(package_status),
    timestamp = format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%S%z"
    )
  )
}

#* Validate an uploaded genomic dataset
#* @param dataset:file Uploaded CSV, TSV, or TXT file
#* @post /validate
function(dataset) {
  tryCatch(
    {
      if (is.null(dataset$datapath)) {
        stop("No dataset was uploaded.")
      }

      extension <- tolower(
        tools::file_ext(dataset$name)
      )

      data <- switch(
        extension,
        csv = readr::read_csv(
          dataset$datapath,
          show_col_types = FALSE,
          name_repair = "unique"
        ),
        tsv = readr::read_tsv(
          dataset$datapath,
          show_col_types = FALSE,
          name_repair = "unique"
        ),
        txt = readr::read_delim(
          dataset$datapath,
          delim = "\t",
          show_col_types = FALSE,
          name_repair = "unique"
        ),
        stop(
          "Supported file formats are CSV, TSV, and TXT."
        )
      )

      gene_input <- tryCatch(
        extract_gene_list(data),
        error = function(error) NULL
      )

      expression_input <- tryCatch(
        prepare_expression_matrix(data),
        error = function(error) NULL
      )

      gene_compatible <- !is.null(gene_input) &&
        length(gene_input$genes) > 0

      gsva_compatible <- !is.null(expression_input) &&
        nrow(expression_input) > 1 &&
        ncol(expression_input) > 1

      list(
        success = TRUE,
        filename = dataset$name,
        rows = nrow(data),
        columns = ncol(data),
        column_names = names(data),
        detected_gene_column = if (gene_compatible) {
          gene_input$column
        } else {
          NULL
        },
        detected_gene_count = if (gene_compatible) {
          length(gene_input$genes)
        } else {
          0
        },
        compatibility = list(
          go = gene_compatible,
          kegg = gene_compatible,
          chea = gene_compatible,
          gsva = gsva_compatible
        )
      )
    },
    error = function(error) {
      list(
        success = FALSE,
        error = conditionMessage(error)
      )
    }
  )
}

#* Run selected real analysis agents
#* @param dataset:file Uploaded CSV, TSV, or TXT file
#* @param agents Comma-separated agents: go,kegg,gsva,chea
#* @param pvalue_cutoff Statistical cutoff
#* @post /analyze
function(
  dataset,
  agents = "go,kegg,gsva,chea",
  pvalue_cutoff = 0.05
) {
  tryCatch(
    {
      if (is.null(dataset$datapath)) {
        stop("No dataset was uploaded.")
      }

      extension <- tolower(
        tools::file_ext(dataset$name)
      )

      data <- switch(
        extension,
        csv = readr::read_csv(
          dataset$datapath,
          show_col_types = FALSE,
          name_repair = "unique"
        ),
        tsv = readr::read_tsv(
          dataset$datapath,
          show_col_types = FALSE,
          name_repair = "unique"
        ),
        txt = readr::read_delim(
          dataset$datapath,
          delim = "\t",
          show_col_types = FALSE,
          name_repair = "unique"
        ),
        stop(
          "Supported file formats are CSV, TSV, and TXT."
        )
      )

      selected_agents <- unique(
        trimws(
          tolower(
            strsplit(
              agents,
              ",",
              fixed = TRUE
            )[[1]]
          )
        )
      )

      allowed_agents <- c(
        "go",
        "kegg",
        "gsva",
        "chea"
      )

      selected_agents <- intersect(
        selected_agents,
        allowed_agents
      )

      if (length(selected_agents) == 0) {
        stop(
          "Select at least one valid agent: go, kegg, gsva, chea."
        )
      }

      pvalue_cutoff <- suppressWarnings(
        as.numeric(pvalue_cutoff)
      )

      if (
        is.na(pvalue_cutoff) ||
        pvalue_cutoff <= 0 ||
        pvalue_cutoff > 1
      ) {
        pvalue_cutoff <- 0.05
      }

      started_at <- Sys.time()

      agent_results <- lapply(
        selected_agents,
        function(agent_name) {
          result <- run_selected_real_agent(
            agent = agent_name,
            data = data,
            pvalue_cutoff = pvalue_cutoff
          )

          list(
            agent = agent_name,
            success = isTRUE(result$success),
            rows = result$rows,
            message = result$message,
            csv = result$csv,
            plot = result$plot
          )
        }
      )

      names(agent_results) <- selected_agents

      finished_at <- Sys.time()

      list(
        success = all(
          vapply(
            agent_results,
            function(result) {
              isTRUE(result$success)
            },
            logical(1)
          )
        ),
        filename = dataset$name,
        rows = nrow(data),
        columns = ncol(data),
        selected_agents = selected_agents,
        pvalue_cutoff = pvalue_cutoff,
        runtime_seconds = as.numeric(
          difftime(
            finished_at,
            started_at,
            units = "secs"
          )
        ),
        results = agent_results
      )
    },
    error = function(error) {
      list(
        success = FALSE,
        error = conditionMessage(error)
      )
    }
  )
}
