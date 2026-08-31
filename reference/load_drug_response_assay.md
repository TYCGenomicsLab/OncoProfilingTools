# Load a PRISM drug response assay

Loads a PRISM primary repurposing response matrix into a
\`DrugResponseAssay\`.

## Usage

``` r
load_drug_response_assay(
  data_path,
  metadata_path = NULL,
  assay_name = "DrugResponse",
  unit = "logfoldchange",
  identifier = NULL,
  ...
)
```

## Arguments

- assay_name:

  A character string specifying the assay name.

- unit:

  A character string specifying the unit of measurement.

- path:

  A character string specifying the file path to the drug response assay
  data.

- compound_metadata_path:

  Optional path to PRISM compound metadata.

- cell_line_metadata_path:

  Optional path to PRISM cell line metadata.

## Value

A \`DrugResponseAssay\`, which inherits from \`SummarizedExperiment\`.

## Details

The PRISM response matrix is expected to contain treatments/compounds as
rows and DepMap model IDs as columns. This already matches the
Bioconductor convention used by \`SummarizedExperiment\`: rows are
features and columns are samples/models.

Compound metadata from \`Extended_Primary_Compound_List.csv\` is stored
in \`rowData()\`. Cell line metadata from \`Cell_Line_Meta_Data.csv\` is
stored in \`colData()\`.
