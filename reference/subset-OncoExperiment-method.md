# Subset an OncoExperiment

Subset an \`OncoExperiment\` by evaluating the expression against the
experiment-level \`colData()\`. This keeps the assay columns and
metadata in sync while filtering models.

## Usage

``` r
# S4 method for class 'OncoExperiment'
subset(object, subset, select = TRUE, ...)
```

## Arguments

- object:

  An \`OncoExperiment\` object to subset.

- subset:

  A logical expression evaluated in the context of \`colData(object)\`.

- select:

  Currently unused. Included for compatibility with the base
  \`subset()\` generic.

## Value

A new \`OncoExperiment\` object containing only the selected models.
