# Create an empty OncoExperiment object

Creates an empty oncology-focused \`MultiAssayExperiment\` subclass.

## Usage

``` r
OncoExperiment(project_name = "OncoExperiment Project")
```

## Arguments

- project_name:

  A character string representing the project name.

## Value

An empty object of class \`OncoExperiment\`.

## Details

Assays and metadata should be added after construction using package
helper functions such as \`load_assays()\` and \`add_assay()\`.

## Examples

``` r
experiment <- OncoExperiment()
```
