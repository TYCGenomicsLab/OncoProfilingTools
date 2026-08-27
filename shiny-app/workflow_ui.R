# Route-level UI components. Scientific execution remains in app.R and the
# existing agent backends; this file only defines the shorter user journey.

biomarker_agent_keys <- setdiff(names(module_meta), "drug")
drug_agent_keys <- "drug"

landing_page_ui <- function() {
  section(
    class = "landing",
    div(
      class = "landing-copy",
      span(class = "eyebrow", "LOCAL · PRIVATE · REPRODUCIBLE"),
      h2("From molecular data to a defensible biological story"),
      p("Choose a focused workflow. Each route combines reproducible analysis, interactive evidence exploration, and a detailed grounded interpretation from Local, Premium, or Compare mode."),
      div(
        class = "landing-highlight-row",
        span("10 scientific agents"), span("Local + Premium AI"), span("Self-contained reports")
      )
    ),
    div(
      class = "route-grid",
      tags$button(
        id = "nav_biomarker", type = "button", class = "btn action-button route-card route-biomarker", `data-val` = "0",
        span(class = "route-icon", "01"),
        span(class = "route-copy", strong("Biomarker Discovery"), tags$small("Gene lists, DEG tables and expression matrices")),
        span(class = "route-arrow", "→")
      ),
      tags$button(
        id = "nav_drug", type = "button", class = "btn action-button route-card route-drug", `data-val` = "0",
        span(class = "route-icon", "02"),
        span(class = "route-copy", strong("Drug Sensitivity"), tags$small("PRISM-style and long-format drug-response tables")),
        span(class = "route-arrow", "→")
      )
    ),
    div(
      class = "landing-principles",
      div(strong("Observed results"), span("Deterministic outputs stay separate from interpretation")),
      div(strong("Provider choice"), span("Ollama stays local; OpenAI requires consent")),
      div(strong("Research use"), span("Reports carry provenance and limitations"))
    )
  )
}

workflow_page_ui <- function(mode) {
  biomarker <- identical(mode, "biomarker")
  title <- if (biomarker) "Biomarker Discovery" else "Drug Sensitivity"

  tagList(
    section(
      class = paste("workflow-hero", paste0("workflow-", mode)),
      div(
        class = "workflow-navigation",
        actionLink("back_home", "← Back to main menu", class = "back-link back-link-primary"),
        span("Home"), span("/"), strong(title)
      ),
      span(class = "eyebrow", if (biomarker) "9 SCIENTIFIC AGENTS" else "DEDICATED PHARMACOGENOMICS ROUTE"),
      h2(title),
      p(if (biomarker) {
        "Run compatible enrichment, regulatory, network, pathway-activity and immune-composition agents."
      } else {
        "Rank observed assay responses independently from biomarker enrichment; response rank alone does not establish mechanism or clinical efficacy."
      }),
      div(
        class = "workflow-capabilities",
        div(strong("01"), span("Validate input")),
        div(strong("02"), span("Run compatible agents")),
        div(strong("03"), span("Explore interactive evidence")),
        div(strong("04"), span("Review integrated interpretation"))
      )
    ),
    div(
      class = "workflow-grid",
      section(
        class = "panel input-panel",
        div(class = "section-heading", span(class = "step", "01"), div(h2("Upload and validate"), p(if (biomarker) "DEG list or expression matrix" else "Drug-response table"))),
        div(
          class = "upload-zone",
          div(class = "upload-icon", if (biomarker) "DEG" else "AUC"),
          div(h3("Choose analysis file"), p("CSV, TSV or TXT · up to 1 GB")),
          fileInput("dataset", label = NULL, accept = c(".csv", ".tsv", ".txt"), buttonLabel = "Choose file", placeholder = "No file selected")
        ),
        uiOutput("dataset_status"),
        uiOutput("dataset_metrics"),
        uiOutput("gene_selection_note")
      ),
      section(
        class = "panel module-panel",
        div(class = "section-heading", span(class = "step", "02"), div(h2("Configure and run"), p("Only compatible agents can be selected"))),
        if (biomarker) uiOutput("biomarker_configuration") else div(
          class = "scientific-note",
          strong("Response direction"),
          "The agent reports the observed assay metric and only describes lower values as greater measured sensitivity when that direction is explicit."
        ),
        if (biomarker) {
          tagList(
            module_selector_group("gene", agent_groups$gene),
            module_selector_group("expression", agent_groups$expression)
          )
        } else {
          module_selector_group("drug", agent_groups$drug)
        },
        tags$details(
          class = "ollama-settings",
          tags$summary("AI interpretation settings · Local, Premium, or Compare"),
          div(
            class = "ollama-settings-body",
            radioButtons(
              "interpretation_provider",
              "Interpretation mode",
              choices = c(
                "Local Ollama · private" = "ollama",
                "OpenAI Premium · external" = "openai",
                "Compare both · side by side" = "compare"
              ),
              selected = "ollama",
              inline = TRUE
            ),
            conditionalPanel(
              "input.interpretation_provider == 'ollama' || input.interpretation_provider == 'compare'",
              div(
                class = "ollama-settings-grid",
                textInput("ollama_host", "Ollama host", value = "http://127.0.0.1:11434"),
                textInput("ollama_model", "Local model", value = "llama3.1:8b"),
                numericInput("ollama_timeout", "Local timeout (seconds)", value = 300, min = 30, max = 600, step = 30)
              ),
              p(class = "provider-privacy-note provider-local-note", "Local mode keeps the structured prompt on this computer and accepts loopback Ollama hosts only.")
            ),
            conditionalPanel(
              "input.interpretation_provider == 'openai' || input.interpretation_provider == 'compare'",
              div(
                class = "openai-premium-panel",
                div(class = "premium-provider-heading", span("PREMIUM PROVIDER"), strong("OpenAI Responses API"), uiOutput("openai_key_status")),
                div(
                  class = "ollama-settings-grid",
                  selectInput(
                    "openai_model",
                    "OpenAI model",
                    choices = c(
                      "GPT-5.6 Sol · maximum quality" = "gpt-5.6-sol",
                      "GPT-5.6 Terra · balanced" = "gpt-5.6-terra",
                      "GPT-5.6 Luna · economical" = "gpt-5.6-luna"
                    ),
                    selected = "gpt-5.6-terra"
                  ),
                  selectInput("openai_reasoning", "Reasoning effort", choices = c("Low" = "low", "Medium" = "medium", "High" = "high"), selected = "medium"),
                  numericInput("openai_timeout", "Premium timeout (seconds)", value = 240, min = 30, max = 600, step = 30)
                ),
                checkboxInput(
                  "openai_data_consent",
                  "I approve sending the structured result digest to the OpenAI API. No raw uploaded file or patient identifier is sent.",
                  value = FALSE
                ),
                p(class = "provider-privacy-note provider-external-note", "The API key is read only from OPENAI_API_KEY on the R server. It is never entered in this page, cached, logged, or included in reports. Requests use store=false.")
              )
            ),
            p("Every mode uses the same deterministic observations, versioned JSON contract, grounding checks, sanitization, terminal-state handling, and reproducible report provenance.")
          )
        ),
        actionButton("run_analysis", paste("Run", title), class = "primary-button", icon = icon("play"))
      )
    ),
    uiOutput("run_status_bar"),
    uiOutput("results_center")
  )
}

build_app_ui <- function() {
  fluidPage(
    tags$head(
      tags$title("OncoProfilingTools"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$link(rel = "stylesheet", href = "styles.css?v=openai-provider-2"),
      tags$link(rel = "stylesheet", href = "pastel.css?v=openai-provider-2"),
      tags$script(src = "status.js?v=compact-progress-1")
    ),
    div(
      class = "app-shell",
      header(
        class = "site-header",
        div(class = "brand-mark", "OP"),
        div(class = "brand-copy", h1("OncoProfilingTools"), p("Local biomarker and pharmacogenomic research")),
        uiOutput("ai_privacy_status", inline = TRUE)
      ),
      main(class = "dashboard", uiOutput("app_content")),
      footer(class = "site-footer", span("OncoProfilingTools · Research use only"), span("Findings require expert and experimental validation"))
    )
  )
}
