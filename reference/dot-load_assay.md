# Load a single assay using a registered loader

Load a single assay using a registered loader

## Usage

``` r
.load_assay(
  assay_type,
  OncoAssay_class,
  loader,
  data_path,
  metadata_path = NULL,
  identifier = NULL,
  name = NULL,
  ...
)
```

## Arguments

- assay_type:

  Character scalar. Project-level assay type.

- OncoAssay_class:

  Character scalar. Expected returned assay class.

- loader:

  Character scalar. Registered loader function name.

- ...:

  Additional arguments passed to the registered loader.

- path:

  Character scalar. Path to the assay file.

- filename:

  Character scalar. File name being loaded.

## Value

An assay object of the expected class.
