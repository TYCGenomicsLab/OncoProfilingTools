# Add an assay experiment to an OncoExperiment

Adds a named assay experiment to the \`experiments()\` of an
\`OncoExperiment\`.

## Usage

``` r
# S4 method for class 'OncoExperiment,SummarizedExperiment'
add_assay(object, assay, name = NULL, overwrite = FALSE, ...)
```

## Arguments

- object:

  An \`OncoExperiment\` object.

- assay:

  An assay object, such as a \`SummarizedExperiment\`,
  \`ExpressionAssay\`, or \`SingleCellExperiment\`.

- name:

  A single character string specifying the experiment name.

- overwrite:

  Logical value indicating whether to overwrite an existing experiment
  with the same name.

- ...:

  Additional arguments. Currently unused.

## Value

The updated \`OncoExperiment\` object.
