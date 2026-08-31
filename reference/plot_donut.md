# Create a donut plot

Creates a donut plot from assay-level metadata stored in an
\`OncoExperiment\`.

## Usage

``` r
plot_donut(object, variable, assay, variable_expr = substitute(variable), ...)
```

## Arguments

- object:

  An \`OncoExperiment\` object.

- variable:

  Metadata field to plot. Can be unquoted, a character string, or an
  enum value such as \`ModelMetadataFields\$OncotreeLineage\`.

- assay:

  A single character string specifying which assay experiment to use.
  For example, \`"Expression"\` or \`"PRISM"\`.

- variable_expr:

  Captured unevaluated variable expression. Internal.

- ...:

  Additional arguments passed to methods.

## Value

A \`ggplot\` object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_donut(exp, ModelMetadataFields$OncotreeLineage, assay = "Expression")
plot_donut(exp, "Sex", assay = "Expression")
plot_donut(exp, OncotreeLineage, assay = "Expression")
} # }
```
