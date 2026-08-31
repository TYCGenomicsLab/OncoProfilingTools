# An enumeration of supported assay types

Helpful for intelligent code completion and avoiding typos in assay
types.

## Usage

``` r
AssayTypes
```

## Examples

``` r
AssayTypes # you can see all supported assays available
#> # A generic enum: 2 members
#>  chr RNA          : RNA
#>  chr DrugResponse : DrugResponse
AssayTypes$RNA # it's just a string at the end of the day
#> [1] "RNA"

```
