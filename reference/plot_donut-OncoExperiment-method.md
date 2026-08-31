# Create a donut plot from model metadata

Generates a donut plot that shows the distribution of categorical
metadata fields from a selected assay in an \`OncoExperiment\`.

## Usage

``` r
# S4 method for class 'OncoExperiment'
plot_donut(
  object,
  variable,
  assay,
  variable_expr = substitute(variable),
  show_na = FALSE,
  ...
)
```

## Arguments

- object:

  An \`OncoExperiment\` object.

- variable:

  Metadata field to plot.

- assay:

  Assay name to use.

- variable_expr:

  Captured unevaluated variable expression. Internal.

- show_na:

  Logical. Whether missing values should be shown as \`"Missing"\`.

- ...:

  Additional arguments passed to \`.apply_onco_plot_theme()\`.

## Value

A \`ggplot\` object.
