# Load assay data into an OncoExperiment object The function will coerce the assay into the corresponding OncoAssay class before being added to the OncoExperiment object.

Load assay data into an OncoExperiment object The function will coerce
the assay into the corresponding OncoAssay class before being added to
the OncoExperiment object.

## Usage

``` r
load_assays(object, assay_type = NULL, data = NULL, metadata = NULL, ...)
```

## Arguments

- object:

  An OncoExperiment object to which the assay will be added.

- assay_type:

  A character string of the type of assay being loaded. Use the
  \`AssayTypes\` enum to easily select a valid assay type.

- data:

  A character string specifying specific file path where the assay data
  is located.

- metadata:

  A character string specifying the path to the metadata file. This is
  optional but highly recommended.

## Value

An updated OncoExperiment with loaded assays in the @assays slot.
