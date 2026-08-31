# OncoExperiment S4 class definition

\`OncoExperiment\` is a \`MultiAssayExperiment\` subclass for oncology
profiling projects.

## Value

An object of class \`OncoExperiment\`.

## Details

It stores multiple assay datasets, such as bulk expression, PRISM drug
response, mutation data, copy number data, and single-cell experiments,
using the standard Bioconductor \`MultiAssayExperiment\` infrastructure.

Additional oncology-specific results, such as fitted models and graphs,
are stored in dedicated slots.

## Slots

- `project_name`:

  A character string representing the project name.

- `version`:

  A character string representing the package/object version.

- `models`:

  A list of fitted models, predictions, or other analytical results.

- `graphs`:

  A list of plots or graph objects related to the experiment.
