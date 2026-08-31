# Create a ProteinExpressionAssay object

Creates a \`ProteinExpressionAssay\`, a \`SummarizedExperiment\`
subclass for protein-coding gene expression data.

## Usage

``` r
RNAAssay(
  data,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  assay_name = "RNA",
  unit = "log2(TPM+1)",
  normalized = TRUE,
  feature_type = "protein_coding_gene",
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

  Name of the primary assay. Defaults to \`"protein_expression"\`.

- unit:

  Character string describing the expression unit, such as
  \`"log2(TPM+1)"\`.

- normalized:

  Logical scalar indicating whether the expression data are normalized.

- feature_type:

  Character string describing the feature type, such as
  \`"protein_coding_gene"\`.

- source:

  Optional character string describing the data source.

- source_file:

  Optional character string storing the source file path.

## Value

A \`ProteinExpressionAssay\` object.
