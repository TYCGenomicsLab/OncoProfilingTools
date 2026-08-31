# Basic Example

## Basic Example \## Installation This package is still under development, so we’ll install it locally from source. I recommend the R `pak` package. It’s becoming the default, is very nice, and handles dependencies well. Then, you can install the `OncoProfilingTools` package from the local source directory like this:

``` r

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pkg_install(".", ask = FALSE)

library(OncoProfilingTools)
```
