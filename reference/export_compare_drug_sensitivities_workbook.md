# Export a drug sensitivity comparison to an Excel workbook

Internal workbook writer used by the \`export()\` S4 method for
\`DrugSensitivityComparison\` objects.

## Usage

``` r
export_compare_drug_sensitivities_workbook(
  comparison_result,
  group1_label,
  group2_label,
  output_file = "drug_sensitivities_comparison.xlsx",
  overwrite = TRUE,
  include_volcano = TRUE,
  n_label_volcano = 10L
)
```

## Arguments

- comparison_result:

  A \`DrugSensitivityComparison\` object.

- group1_label:

  Character scalar used for display names in the workbook.

- group2_label:

  Character scalar used for display names in the workbook.

- output_file:

  Character scalar path to the \`.xlsx\` file to create.

- overwrite:

  Logical value; overwrite an existing file if \`TRUE\`.

- include_volcano:

  Logical value; if \`TRUE\`, add a volcano plot sheet.

- n_label_volcano:

  Integer scalar controlling how many features to label.

## Value

Invisibly returns \`output_file\`.
