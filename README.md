<p align="center" width="100%">
    <img style="max-width: 200px; width: 50%; height: auto;" src="https://avatars.githubusercontent.com/u/291118088?s=200&amp;v=4" alt="TYC Genomics Lab logo"/>
    <img style="max-width: 200px; width: 50%; height: auto;" src="https://images.seeklogo.com/logo-png/37/2/virginia-commonwealth-university-vcu-logo-png_seeklogo-376239.png" alt="Virginia Commonwealth University logo"/>
</p>

---

# OncoProfilingTools

This package is currently under active development and is not ready for production or clinical use.

[![pkgdown.yaml](https://github.com/TYCGenomicsLab/OncoProfilingTools/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/TYCGenomicsLab/OncoProfilingTools/actions/workflows/pkgdown.yaml)

## Shiny research interface

The development branch also contains a local Shiny interface for biomarker enrichment and drug-sensitivity workflows. Its Results Center separates deterministic observations from optional model interpretation and produces a compact combined HTML report.

Interpretation modes are:

- **Rules-based** — deterministic summaries; no AI service or API key.
- **Local Ollama** — sends the structured result digest to a loopback Ollama server.
- **OpenAI Premium** — asks for an API key only in the local browser after this mode is selected.
- **Compare both** — shows Ollama and OpenAI interpretations from the same result digest.

The terminal never asks for an OpenAI key. A key entered in the browser is held only for the current Shiny session and is not written to reports, logs, caches, or settings files.

Start the application:

```bash
Rscript --vanilla -e 'shiny::runApp("shiny-app", host="127.0.0.1", port=3838, launch.browser=TRUE)'
```

Run the complete reviewer verification suite:

```bash
Rscript --vanilla scripts/verify_pr.R
```

See [TESTING.md](TESTING.md) for the click-by-click testing workflow and [RESULTS_CENTER_ARCHITECTURE.md](RESULTS_CENTER_ARCHITECTURE.md) for the interpretation contract, prompt locations, safety boundaries, and extension points.

> Research use only. Results require independent statistical, biological, and clinical review.
