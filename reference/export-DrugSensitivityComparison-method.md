# Export a \`DrugSensitivityComparison\` object

Export a \`DrugSensitivityComparison\` object

## Usage

``` r
# S4 method for class 'DrugSensitivityComparison'
export(
  object,
  output_file = "drug_sensitivities_comparison.xlsx",
  overwrite = TRUE,
  include_volcano = TRUE,
  n_label_volcano = 10L,
  group1_label = object@group1_label %||% "Group 1",
  group2_label = object@group2_label %||% "Group 2",
  ...
)
```

## Arguments

- object:

  A \`DrugSensitivityComparison\` object.

- output_file:

  Character scalar path to the \`.xlsx\` file to create.

- overwrite:

  Logical value; overwrite an existing file if \`TRUE\`.

- include_volcano:

  Logical value; if \`TRUE\`, add a volcano plot sheet.

- n_label_volcano:

  Integer scalar controlling how many features to label.

- group1_label:

  Optional override for the workbook label of group 1.

- group2_label:

  Optional override for the workbook label of group 2.

- ...:

  Additional arguments passed through to the workbook exporter.

## Value

Invisibly returns \`output_file\`.
