# Create an ExpressionAssay object

Creates an \`ExpressionAssay\`, a \`SummarizedExperiment\` subclass for
expression data.

## Usage

``` r
ExpressionAssay(
  data,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  assay_name = "expression",
  unit = "log2(TPM+1)",
  normalized = TRUE,
  feature_type = "transcript",
  source = NA_character_,
  source_file = NA_character_
)
```

## Arguments

- data:

  A numeric matrix-like object with features as rows and samples/models
  as columns.

- rowData:

  Optional feature-level metadata. Must have one row per feature.

- colData:

  Optional sample/model-level metadata. Must have one row per
  sample/model.

- metadata:

  Optional assay-level metadata list.

- assay_name:

  Name of the primary assay. Defaults to \`"expression"\`.

- unit:

  Character string describing the expression unit, such as
  \`"log2(TPM+1)"\`.

- normalized:

  Logical scalar indicating whether the expression data are normalized.

- feature_type:

  Character string describing the feature type, such as \`"gene"\` or
  \`"transcript"\`.

- source:

  Optional character string describing the data source.

- source_file:

  Optional character string storing the source file path.

## Value

An \`ExpressionAssay\` object.
