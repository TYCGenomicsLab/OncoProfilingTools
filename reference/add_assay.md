# Add an assay to an object

Adds an assay object to a container object.

## Usage

``` r
add_assay(object, assay, name = NULL, overwrite = FALSE, ...)
```

## Arguments

- object:

  The object to which the assay will be added.

- assay:

  The assay object to add.

- name:

  Optional character string specifying the assay name. If \`NULL\`, the
  assay object's existing \`name\` slot is used.

- overwrite:

  Logical value indicating whether to overwrite an existing assay with
  the same name.

- ...:

  Additional arguments passed to methods.

## Value

The updated object.
