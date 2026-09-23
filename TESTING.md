# Reviewer testing guide

This guide covers the Shiny frontend, the ten analysis modules, the interpretation providers, and the corrected CMS4 inputs supplied for testing.

## 1. Start the application

From the repository root in the VS Code terminal, run:

```bash
Rscript --vanilla -e 'shiny::runApp("shiny-app", host="127.0.0.1", port=3838, launch.browser=TRUE)'
```

The app should open in Chrome at `http://127.0.0.1:3838`. This command must not request an OpenAI API key.

For Local Ollama or Compare mode, start Ollama separately if it is not already running:

```bash
ollama list
ollama pull llama3.1:8b
ollama serve
```

On macOS, `ollama serve` may report that port `11434` is already in use because the Ollama desktop application has already started the service. In that case, do not start a second service.

## 2. Test the frontend

1. Confirm the landing page has separate **Biomarker Discovery** and **Drug Sensitivity** routes.
2. Open **Biomarker Discovery** and upload a CSV, TSV, or TXT file.
3. Confirm the upload summary reports the filename, dimensions, header detection, gene column, identifier mapping, and automatic gene-selection rule.
4. Confirm compatible module cards are enabled and checked automatically, while incompatible cards remain disabled.
5. Confirm the interpretation selector includes:
   - Rules-based · no AI
   - Local Ollama · private
   - OpenAI Premium · external
   - Compare both · side by side
6. Select Rules-based or Local Ollama. The OpenAI API-key field must not be visible.
7. Select OpenAI Premium or Compare both. The session-only password field and consent checkbox must appear in the browser.
8. Run compatible agents and confirm the compact status strip reaches a terminal completed or failed state for every selected agent.
9. Confirm the Results Center provides per-agent evidence, interpretation, CSV/report downloads, a Combined HTML Report, and a complete ZIP bundle.
10. Open the Combined HTML Report and confirm repetitive technical tables are collapsed or omitted, while exact prompts/rules and full-precision downloadable artifacts remain available.

## 3. Test the corrected CMS4 datasets

Use the corrected files containing both `entrez_gene_id` and `gene_symbol`.

| Input | Expected automatic selection | Expected compatible modules |
|---|---:|---|
| `CCLE_CRC_volcano_right_higher_in_CMS4.csv` | 586 higher-in-CMS4 genes | GO, KEGG, Reactome, WikiPathways, STRING, Hallmark, ChEA |
| `CCLE_CRC_volcano_left_lower_in_CMS4.csv` | 489 lower-in-CMS4 genes | GO, KEGG, Reactome, WikiPathways, STRING, Hallmark, ChEA |
| `CCLE_CRC_all_genes_ranked_CMS4_vs_Other.csv` | 1,075 genes using `FDR <= 0.05` and `abs(log2FC_CMS4_vs_Other) >= 1` | Seven gene-list modules only |
| `CCLE_CRC_gene_limma_CMS4_vs_Other_full_results.csv` | 1,075 genes using `adj.P.Val <= 0.05` and `abs(logFC) >= 1` | Seven gene-list modules only |

For all four files, GSVA and Immune Deconvolution must remain disabled because differential-expression statistics are not sample-level expression columns. Drug Sensitivity must also remain disabled.

Run the higher- and lower-in-CMS4 files separately so the biological direction is preserved. The full ranked file is currently converted into the same thresholded gene-list input used for over-representation analysis; it is not a true preranked GSEA implementation.

## 4. Test each interpretation mode

### Rules-based

Choose **Rules-based · no AI** and run the selected modules. No model service or API key is required. The report should explicitly state that no language-model prompt ran.

### Local Ollama

Choose **Local Ollama · private**, keep the loopback host and installed model, and run. Only the structured result digest should be sent to the local Ollama service.

### OpenAI Premium

Choose **OpenAI Premium · external**. The password field appears only now. Paste the OpenAI Platform key into the browser field, select the consent checkbox, and run.

Do not enter or export the key in the VS Code terminal. The key is used only by the current local Shiny session, is passed directly to the background interpretation process, and is not written to settings files, reports, logs, caches, or Git.

### Compare both

Start Ollama, choose **Compare both · side by side**, paste the OpenAI key in the browser field, approve consent, and run. The report should show both providers and state that they received the same structured digest and prompt contract.

## 5. Module-specific expectations

| Module | Test input | Expected result |
|---|---|---|
| GO | CMS4 gene list | Biological Process enrichment CSV and dot plot |
| KEGG | CMS4 gene list | Human KEGG pathway CSV and dot plot |
| Reactome | CMS4 gene list | Reactome pathway CSV and plot |
| WikiPathways | CMS4 gene list | WikiPathways enrichment CSV and plot |
| STRING | CMS4 gene list | Hub table, interaction edge list, static network, and connected 3D report view |
| Hallmark | CMS4 gene list | Hallmark enrichment CSV and plot |
| ChEA | CMS4 gene list | ChEA 2022 regulator table and plot |
| GSVA | Real expression matrix with genes plus at least two sample columns | Hallmark score matrix and heatmap |
| Immune Deconvolution | Real expression matrix with genes plus at least two sample columns | Immune-composition table and heatmap |
| Drug Sensitivity | Long-form drug/response table or wide PRISM-style table | Ranked compound table and response plot |

KEGG, WikiPathways, STRING, and ChEA can contact external reference services. A network or upstream-service failure should be reported as a failed module without stopping successful modules or fabricating results.

## 6. Run automated verification all at once

From the repository root, run the complete verification with one command:

```bash
Rscript --vanilla scripts/verify_pr.R
```

This runs the test suite, parses every R source file, checks the frontend JavaScript, and validates the Git diff. Provider-contract tests use mocks and do not spend OpenAI credits. Live analysis modules are tested manually through the browser with the appropriate input type.

## 7. Troubleshooting

- **Port already in use:** stop the previous Shiny process or change `port=3838` to another port such as `3841`.
- **Ollama unavailable:** run `ollama list`; start the desktop app or local service if needed.
- **OpenAI key rejected:** replace the value in the browser password field and run again. Do not restart the terminal or export the key.
- **A gene-list module is disabled:** confirm the upload contains a recognized gene-symbol, Entrez-ID, or Ensembl-ID column and at least two mapped genes.
- **GSVA or Immune is disabled:** use a true sample-by-gene or gene-by-sample expression matrix, not a DEG statistics table.
