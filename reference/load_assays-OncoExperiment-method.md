# Load assay data into an OncoExperiment object

This function loads and parses the specified data and metadata files and
appends to the @assays slot of the \`OncoExperiment\` object. You can
view available assay types by printing \`AssayTypes\`.

## Usage

``` r
# S4 method for class 'OncoExperiment'
load_assays(
  object,
  assay_type = NULL,
  data = NULL,
  metadata = NULL,
  identifier = NULL,
  name = NULL,
  overwrite = FALSE,
  ...
)
```

## Arguments

- object:

  An \`OncoExperiment\` object to which the assay will be added.

- assay_type:

  A character string of the type of assay being loaded. Use the
  \`AssayTypes\` enum to easily select a valid assay type.

- data:

  A character string specifying the path to the data file. This is
  required.

- metadata:

  A character string specifying the path to the metadata file. This is
  optional but highly recommended.

- identifier:

  A character string of the column name where the unique identifier is
  found that links your data to the metadata. If not supplied, the first
  column is used.

- name:

  A character string of a unique name for this assay. If not supplied
  and an assay of the same type exists, they will be named sequentially.

- overwrite:

  A logical value indicating whether to overwrite an existing assay with
  the same name.

## Value

The updated \`OncoExperiment\` object.
