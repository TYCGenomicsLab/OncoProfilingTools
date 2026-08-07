# Results Center architecture

## Scope

The Results Center leaves the ten scientific workers and their canonical CSV/plot outputs intact. A separate adapter reads those outputs and supplies presentation, local interpretation, and synthesis. This keeps analysis execution independent from model availability.

## Layout and grouping

- Gene-list agents: GO, KEGG, Reactome, WikiPathways, STRING, Hallmark, and ChEA.
- Expression profiling: GSVA and Immune Deconvolution.
- Drug Sensitivity: a dedicated section.
- Every agent tab uses a balanced two-column desktop grid. Visualization and result table are stacked in the left evidence column; a sticky, independently scrollable biological interpretation occupies the right column. The grid collapses to one column on narrower screens.

Ordinary numeric scores are formatted in the browser to two decimal places. Probability columns such as p-values, adjusted p-values, q-values, and FDR use scientific notation so small values never appear as a misleading `0.00`. Download handlers copy the canonical CSV rather than serializing the display table, so source precision is retained.

## Gene-set preparation and upstream validation

DEG-style uploads are reduced to an explicit, bounded gene signature before any of the seven gene-list workers run. Adjusted p-value/FDR/q-value columns use a 0.05 cutoff, detected effect columns use an absolute cutoff of 1, and no more than 2,000 strongest rows are retained. The UI discloses the selected and original gene counts. A large table with no ranking statistic is rejected rather than silently using a near-whole-genome background as the query.

ChEA treats the Enrichr response as untrusted upstream data. A valid response must contain `Term`, `Adjusted.P.value`, and `Combined.Score`; expired or malformed payloads are retried up to three times and are never written as result CSVs. STRING writes both its hub table and complete interaction edge list, then renders the strongest connected subnetwork (or a hub-degree fallback when `igraph` is unavailable).

## Standard agent exchange

`build_agent_exchange(agent_id, data)` returns a versioned list with:

- `agent_id` and analysis `domain`;
- `row_count`, `column_count`, and the detected label column;
- up to eight representative findings selected by significance/rank/score when available;
- numeric minimum, median, maximum, and finite-value count calculated across complete columns;
- broad biological programs detected across every result label; and
- a coverage note that states how the complete result set contributed.

The exchange is intentionally concise so agents can share evidence without passing large tables into a prompt. Add future agent adapters by extending `agent_domain()`, `result_label_column()`, and, if needed, `rank_result_rows()`; the synthesis layer does not need direct dependencies on worker functions.

## Drug–pathway hook

`build_drug_pathway_bridge(exchanges)` activates when Drug Sensitivity and at least one pathway/systems agent are present. It publishes ranked compounds, pathway-agent IDs, and detected programs to the synthesis layer. It deliberately makes no drug-mechanism or causal inference. A future matched-sample association service can consume the versioned bridge without changing Drug Sensitivity or GSVA worker code.

## Local Ollama flow

1. Result file signatures change as workers finish.
2. The app waits until the selected run is complete.
3. One exchange is built for every selected agent.
4. The app immediately publishes a deterministic full-result summary and starts a separate local R worker.
5. The worker sends one structured JSON request to `/api/generate` on the configured Ollama host while the Shiny session remains responsive.
6. The response is schema-checked and text is rendered as escaped UI content.
7. Missing, invalid, disabled, or timed-out responses retain the deterministic rule bundle with a concise user-facing message; terminal control codes and raw curl diagnostics are not exposed.
8. The bundle is cached by result signatures, model, host, timeout, response length, and run state.

Only loopback URLs matching `localhost`, `127.0.0.1`, or `[::1]` are accepted. Result contents are marked as untrusted data in the prompt, and model output is never inserted as raw HTML.

Default settings:

```text
Host:  http://127.0.0.1:11434
Model: llama3.1:8b
Timeout: 180 seconds
Maximum response: 1800 tokens
```

Environment overrides are available for non-UI uses:

```bash
export ONCOPROFILING_OLLAMA_HOST=http://127.0.0.1:11434
export ONCOPROFILING_OLLAMA_MODEL=llama3.1:8b
export ONCOPROFILING_OLLAMA_TIMEOUT=180
export ONCOPROFILING_OLLAMA_NUM_PREDICT=1800
```

## Plot refresh

Result polling watches file size and high-resolution modification time. When a plot is rendered, `plot_cache_token()` also includes an MD5 content fingerprint. If `base64enc` is installed, the image is served inline from the current file contents; otherwise the fingerprint is appended to the local resource URL. Both paths prevent an overwritten plot from reusing a stale browser cache entry.

## Verification commands

From the repository root:

```bash
Rscript -e 'testthat::test_dir("tests/testthat", reporter="summary")'
Rscript -e 'files <- c(list.files("R", "[.]R$", full.names=TRUE), list.files("shiny-app", "[.]R$", full.names=TRUE)); invisible(lapply(files, parse)); cat(length(files), "R files parsed\n")'
node --check shiny-app/www/status.js
Rscript -e 'old <- setwd("shiny-app"); on.exit(setwd(old)); app <- source("app.R", local=new.env())$value; stopifnot(inherits(app, "shiny.appobj")); cat("Shiny app source smoke test passed\n")'
```
