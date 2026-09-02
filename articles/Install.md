# Install

## How to Install

### Using `pak`

OncoProflingTools (OPT) is currently available from our [GitHub
repository](https://github.com/oncoproflingtools/oncoproflingtools)
while in active development. In the future, we plan to make it available
on CRAN or Bioconductor.

Ensure [R](https://www.r-project.org/) version 4.0 or greater is
installed. To install OPT, run the following command in the R console:

``` r

install.packages("pak") # alternatively, use native R package manager
pak::pkg_install("TYCGenomicsLab/OncoProfilingTools") # will fetch from GitHub
library(OncoProfilingTools)
```

### Dependencies

It is highly recommended you also install and load
[dplyr](https://dplyr.tidyverse.org/),
[MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment/),
and
[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/)
to make use of the full functionality of OPT. To install these packages,
run:

``` r

# pak will automatically know to use CRAN or Bioconductor
pak::pkg_install(c("dplyr", "MultiAssayExperiment", "SummarizedExperiment"))
library("dplyr")
library("MultiAssayExperiment")
library("SummarizedExperiment")
```

`Pak` will automatically install any dependencies required by OPT. If
you recieve any issues about missing packages, please install them
manually.
