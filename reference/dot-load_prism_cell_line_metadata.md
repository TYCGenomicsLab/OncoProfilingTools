# Load PRISM cell line metadata

Loads and aligns PRISM cell line metadata to drug response assay column
names.

## Usage

``` r
.load_prism_cell_line_metadata(cell_line_metadata_path = NULL, model_ids)
```

## Arguments

- cell_line_metadata_path:

  Path to the PRISM cell line metadata file, or \`NULL\`.

- model_ids:

  Character vector of model IDs, usually assay column names.

## Value

An \`S4Vectors::DataFrame\` aligned to \`model_ids\`.

## Details

The PRISM cell line metadata is expected to contain a \`depmap_id\`
column matching the model IDs in the drug response matrix columns.
Because PRISM metadata may contain multiple rows per \`depmap_id\`
across pools, cultures, and screens, duplicated \`depmap_id\` rows are
collapsed to one row per model.
