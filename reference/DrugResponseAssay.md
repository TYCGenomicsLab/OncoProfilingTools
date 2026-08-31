# Create a DrugResponseAssay object

Creates a \`DrugResponseAssay\`, a \`SummarizedExperiment\` subclass for
drug response data.

## Usage

``` r
DrugResponseAssay(
  data,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  assay_name = "response",
  unit = "logfoldchange",
  source = NA_character_,
  source_file = NA_character_
)
```

## Arguments

- data:

  A numeric matrix-like object with treatments/drugs as rows and
  samples/models as columns.

- rowData:

  Optional treatment/compound-level metadata. Must have one row per
  treatment/drug.

- colData:

  Optional sample/model-level metadata. Must have one row per
  sample/model.

- metadata:

  Optional assay-level metadata list.

- assay_name:

  Name of the primary assay. Defaults to \`"response"\`.

- unit:

  Character string describing the response unit.

- source:

  Optional character string describing the data source.

- source_file:

  Optional character string storing the source file path.

## Value

A \`DrugResponseAssay\` object.
