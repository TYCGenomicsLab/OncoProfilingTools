# Load PRISM compound metadata

Loads and aligns PRISM compound metadata to drug response assay row
names.

## Usage

``` r
.load_prism_compound_metadata(compound_metadata_path = NULL, treatment_ids)
```

## Arguments

- compound_metadata_path:

  Path to the PRISM compound metadata file, or \`NULL\`.

- treatment_ids:

  Character vector of treatment IDs, usually assay row names.

## Value

An \`S4Vectors::DataFrame\` aligned to \`treatment_ids\`.

## Details

The PRISM compound list is expected to contain an \`IDs\` column
matching the treatment IDs in the response matrix. If the compound
metadata contains duplicated \`IDs\`, rows are collapsed to one row per
\`IDs\` by concatenating unique non-missing values with \`"; "\`.
