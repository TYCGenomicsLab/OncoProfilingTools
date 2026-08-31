# Load DepMap Model.csv metadata

Loads and aligns DepMap \`Model.csv\` metadata to assay column names.

## Usage

``` r
.load_model_metadata(model_metadata_path = NULL, model_ids)
```

## Arguments

- model_metadata_path:

  Path to \`Model.csv\`, or \`NULL\`.

- model_ids:

  Character vector of model IDs, usually assay column names.

## Value

An \`S4Vectors::DataFrame\` aligned to \`model_ids\`.
