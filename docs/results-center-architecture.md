# Results Center architecture

## Scope

The opening page routes Biomarker Discovery and Drug Sensitivity independently. The Results Center leaves the ten scientific workers and their canonical CSV/plot outputs intact. A separate adapter reads those outputs and supplies presentation, deterministic observation, local interpretation, reporting, and synthesis. This keeps analysis execution independent from model availability.

## Layout and grouping

- Gene-list agents: GO, KEGG, Reactome, WikiPathways, STRING, Hallmark, and ChEA.
- Expression profiling: GSVA and Immune Deconvolution.
- Drug Sensitivity: a dedicated section.
- Every agent tab uses a balanced two-column desktop grid. A professional interactive horizontal bar view occupies the left evidence column and a sticky, independently scrollable biological interpretation occupies the right. Generic 3D rank profiles and repetitive preview tables are omitted. STRING uniquely adds a mouse-rotatable 3D network because its edges represent retrieved protein associations. The grid collapses to one column on narrower screens.

Ordinary numeric scores are formatted in the browser to two decimal places. Probability columns such as p-values, adjusted p-values, q-values, and FDR use scientific notation so small values never appear as a misleading `0.00`. Download handlers copy the canonical CSV rather than serializing the display table, so source precision is retained.

For quanTIseq immune results, the `Other`/`uncharacterized cell` row is treated as an unresolved residual compartment. It remains in the canonical CSV, but the visualization and interpretation adapter report it separately rather than ranking it as a named immune-cell population. The immune heatmap scales colors across named populations only and states whether all samples or a bounded variable subset are displayed.

## Gene-set preparation and upstream validation

Headerless one-column human Ensembl lists and named HGNC/Ensembl/Entrez columns are detected explicitly. Version suffixes are removed before `org.Hs.eg.db` mapping; unique-input, mapped, unmapped, output-symbol, duplicate-removal, rate, and example metadata travel with the run. DEG-style uploads are reduced to an explicit, bounded gene signature before any of the seven gene-list workers run. IAN-style `NominalSign` and `gene_type` fields expose focused group and protein-coding controls. Adjusted p/FDR and effect thresholds plus the maximum gene count are visible configuration rather than hidden constants. A large table with no ranking statistic is rejected rather than silently using a near-whole-genome background as the query.

ChEA treats the Enrichr response as untrusted upstream data. A valid response must contain `Term`, `Adjusted.P.value`, and `Combined.Score`; expired or malformed payloads are retried up to three times and are never written as result CSVs. STRING writes both its hub table and complete interaction edge list, then renders the strongest connected subnetwork (or a hub-degree fallback when `igraph` is unavailable).

## Standard agent exchange

`build_agent_exchange(agent_id, data)` returns a versioned list with:

- `agent_id` and analysis `domain`;
- `row_count`, `column_count`, and the detected label column;
- up to eight representative findings selected by significance/rank/score when available;
- numeric minimum, median, maximum, and finite-value count calculated across complete columns;
- broad biological programs detected across every result label; and
- a coverage note that states how the complete result set contributed.

The exchange is intentionally concise so agents can share evidence without passing large tables into a prompt. The output contract is `4.3-ian-evidence-integrated-report`. It adapts the IAN combined-prompt sequence into versioned JSON: individual-agent review, cross-agent integration, groundedness checks, candidate regulatory/network interpretation, hypotheses, validation priorities, high-level synthesis, and a concise title. Every agent entry receives `observed_results` and named key findings rebuilt deterministically from its own exchange; the model cannot supply or overwrite these fields. Model summaries, biological context, and hypotheses must mention a ranked observation from that same agent or the computed agent interpretation is retained. Drug mechanisms and targets absent from the observed assay table are rejected. Add future agent adapters by extending `agent_domain()`, `result_label_column()`, and, if needed, `rank_result_rows()`; the synthesis layer does not need direct dependencies on worker functions.

## Drug–pathway hook

`build_drug_pathway_bridge(exchanges)` activates when Drug Sensitivity and at least one pathway/systems agent are present. It publishes ranked compounds, pathway-agent IDs, and detected programs to the synthesis layer. It deliberately makes no drug-mechanism or causal inference. A future matched-sample association service can consume the versioned bridge without changing Drug Sensitivity or GSVA worker code.

## Provider-neutral interpretation flow

1. Result file signatures change as workers finish.
2. The app waits until the selected run is complete.
3. One exchange is built for every selected agent.
4. The app immediately publishes deterministic observed results and starts a separate R worker. While it runs, the UI shows a dedicated interpretation-building state rather than a provisional computed narrative.
5. The worker routes the same structured digest and JSON contract to local Ollama, the OpenAI Responses API, or both providers according to the selected mode while the Shiny session remains responsive.
6. The response must match the exact top-level contract version, is schema-checked and sanitized, and text is rendered as escaped UI content. OpenAI receives a closed strict JSON schema materialized for exactly the selected agent IDs; no dynamic `additionalProperties` map is sent to Structured Outputs.
7. Missing, invalid, disabled, or timed-out responses retain the deterministic computed-summary bundle with a provider-specific terminal unavailable, timed-out, or error state; terminal control codes and raw HTTP diagnostics are not exposed.
8. A parent-process watchdog terminates a worker that outlives its provider-aware deadline, even if an HTTP client or model process hangs.
9. Terminal success and failure bundles are persisted atomically by result signature, provider, models, host, timeouts, response lengths, consent state, and contract version. A matching completed bundle is restored instead of starting duplicate generation.
10. Static HTML exports never embed a loading/generating state. An early export contains the complete deterministic interpretation; an export after completion uses the persisted provider result. Compare mode adds side-by-side provider provenance, timing, tokens, estimated cost, and agent summaries.

## Combined report contract

`production_reporting.R` creates a light, self-contained HTML report with the input filename in the run summary, MD5 identity, header detection, explicit interpretation state/contract, and a responsible-interpretation notice before the scientific narrative. The Methods section integrates DEG selection provenance, experimental-comparison notes, identifier-mapping quality, gene-universe disclosure, and package versions. Report results use one precise interactive horizontal bar chart per agent; repetitive result-preview tables and generic 3D rank views are intentionally omitted because the complete full-precision CSV files ship in the results bundle. STRING uniquely retains a connected 3D interaction network because its retrieved edges are biologically meaningful. A side guide reports considered nodes, eligible edges, top hub degrees, encodings, mouse controls, and the non-causality boundary. The IAN-style synthesis adds a biological/cellular interpretation, deterministic pathway member-gene overlap, ChEA regulator rationale, a visible STRING network with hub rationale, and one consolidated result-grounded hypothesis. Biomarker Discovery and Drug Sensitivity remain distinct, and the artifact manifest is the final report table. Plotly JavaScript is embedded locally so report interaction does not require a CDN. The ZIP bundle adds `artifact_manifest.csv` and `gene_identifier_mapping.csv` when mapping was performed.

Only loopback Ollama URLs matching `localhost`, `127.0.0.1`, or `[::1]` are accepted. OpenAI mode reads `OPENAI_API_KEY` from the R process environment, requires explicit UI consent, sends no original upload or patient identifier, and sets `store=false`. The key is not serialized into worker settings, browser state, caches, logs, or reports. Result contents are marked as untrusted data in the prompt, and model output is never inserted as raw HTML. Failed calls retain only a sanitized provider message, HTTP status, and request ID for troubleshooting. Provider comparison views never label the deterministic computed summary as model-authored output when that provider did not complete.

Default settings:

```text
Host:  http://127.0.0.1:11434
Model: llama3.1:8b
Timeout: 300 seconds
Maximum response: 4096 tokens
```

Environment overrides are available for non-UI uses:

```bash
export ONCOPROFILING_OLLAMA_HOST=http://127.0.0.1:11434
export ONCOPROFILING_OLLAMA_MODEL=llama3.1:8b
export ONCOPROFILING_OLLAMA_TIMEOUT=300
export ONCOPROFILING_OLLAMA_NUM_PREDICT=4096
export OPENAI_API_KEY=your-platform-api-key
export ONCOPROFILING_OPENAI_MODEL=gpt-5.6-terra
export ONCOPROFILING_OPENAI_REASONING=medium
export ONCOPROFILING_OPENAI_TIMEOUT=240
export ONCOPROFILING_OPENAI_MAX_OUTPUT=12000
```

OpenAI Platform API access and billing are separate from the local Ollama path. Do not send identifiable patient data or protected health information without the institutional approvals and provider terms required for that data. Any displayed API cost is an estimate based on rates encoded in the app version; current platform pricing remains authoritative.

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
