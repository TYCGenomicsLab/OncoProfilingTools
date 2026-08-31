# Load a DepMap protein-coding gene expression assay

Loads a DepMap protein-coding gene expression file into a \`RNAAssay\`.

## Usage

``` r
load_RNA_assay(
  data_path,
  metadata_path = NULL,
  assay_name = "Protein Expression",
  unit = "log2(TPM+1)",
  normalized = TRUE,
  id_col = "ModelID",
  feature_type = "protein_coding_gene",
  identifier = NULL,
  ...
)
```

## Arguments

- assay_name:

  A character string specifying the assay name used in metadata.

- unit:

  A character string specifying the expression unit.

- normalized:

  A logical value indicating whether the expression data are normalized.

- id_col:

  A character string specifying the expression file column to use as
  model/sample identifiers.

- feature_type:

  A character string describing the expression feature type.

- path:

  A character string specifying the path to the expression assay file.

- model_metadata_path:

  Optional path to \`Model.csv\`. If \`NULL\`, the loader searches for
  \`Model.csv\` beside \`path\`.

## Value

A \`RNAAssay\`, which inherits from \`SummarizedExperiment\`.

## Details

The DepMap expression file is read, filtered to default model entries,
and converted into a Bioconductor-style expression matrix with features
as rows and models as columns. If \`Model.csv\` is present beside the
assay file, it is loaded and aligned into \`colData()\`.
