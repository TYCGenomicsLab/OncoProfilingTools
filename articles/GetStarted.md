# Getting started

## Understanding the structure

Before starting to interact with the package you should understand the
structure of the objects you’ll be working with. OPT borrows and expands
from well-known Bioconductor packages for dealing with large-scale
experimental data. The following two objects are the most important to
understand:

1.  **`OncoExperiment`** - a class that inherits from
    [MultiAssayExperiment](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html).
    It contains all the information needed to perform downstream
    analyses in the form of assays.
2.  **`*Assay`** - a class that inherits from
    [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html).
    Different types of assays exist to support loading various sources
    with different structures. Some examples include:
    [`RNAAssay`](https://tycgenomicslab.github.io/OncoProfilingTools/articles/),
    [`DrugResponseAssay`](https://tycgenomicslab.github.io/OncoProfilingTools/articles/),
    and
    [`PRISMAssay`](https://tycgenomicslab.github.io/OncoProfilingTools/articles/).

## `OncoExperiment`

The `OncoExperiment` is the core object that contains all the
information. It harmonizes all the assays (or experiments) together and
maps them to the samples. OncoExperiment is primarily a
`MultiAssayExperiment` (`MAE`) with a few additional methods to deal
with oncogenic data. Although this guide will summarize the structure to
a degree, the authors of `MAE` have created some wonderful [quick start
guides](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/QuickStartMultiAssay.html)
and
[cheatsheets](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment_cheatsheet.html)
with in-depth information. Please refer to these if you ever need help
navigating an `OncoExperiment` you’ve created.

> ![](assets/mae.png)  
> Credit: [“Software For The Integration of Multi-Omics Experiments In
> Bioconductor.” Marcel Ramos et.
> al.](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html)

### Slots

A quick overview of the slots you should know about in an
`OncoExperiment`:

- **`project_name`** - The name of the project for your own reference.
- **`ExperimentList`** - Contains all ID-based experimental data. An
  instance of `ExperimentList`. Access with
  `experiments(OncoExperiment)`.
- **`colData`** - Contains all sample-level metadata. An instance of
  `DataFrame`. Access with `colData(OncoExperiment)` or `$`.
- **`sampleMap`** - Helps to relate experiment-specific sample naming or
  replicate observations to the row names in `colData`. Access with
  `sampleMap(OncoExperiment)`.
- **`metadata`** - Storing any additional metadata about the study. This
  is free to include anything. We store some basic information as a
  `list`. Individual assays can store additional metadata in their own
  `metadata` slot. Access with `metadata(OncoExperiment)`.

This quick start guide by the Waldron Lab has a lot more great
information: [Quick Start
Guide](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/QuickStartMultiAssay.html).
There is also a lot of useful examples showing these accessors in use.

### Creating an `OncoExperiment`

``` r

exp <- OncoExperiment(project_name = "Murine RNA-Seq 04")
exp
```

    ## An OncoExperiment object
    ## Project Name: Murine RNA-Seq 04 
    ## Version: 0.1.0 
    ## Experiments: none

### Subsetting

To demonstrate accessing the data in an `OncoExperiment`, let’s
construct a simple example with an `RNAAssay`. We’ve created an example
file with a 50x50 table. Rows contain unique sample names and colums
contain unique gene names. The values are random numbers between 1 and
10. We can peek at it to demonstrate how it is structured:

``` r

head(read.csv("assets/rna_example.csv", header = TRUE, sep = ","))
```

    ##   ModelID TSPAN6..7105. TNMD..64102. DPM1..8813. SCYL3..57147. FIRRM..55732.
    ## 1 VQM-732      4.956577     0.000000    7.577648      3.179411      4.765742
    ## 2 0Xe-933      4.955015     0.617117    7.333933      2.782935      3.735371
    ## 3 dFL-450      3.421952     0.000000    7.546069      2.615880      4.476233
    ## 4 EcW-351      5.196729     0.000000    6.362268      2.144996      3.087183
    ## 5 nkn-474      4.651643     0.000000    5.946408      2.454515      1.852111
    ## 6 yBR-806      4.336705     0.000000    6.879387      2.262824      3.256491

Now let’s load the data into an OncoExperiment:

``` r

example <- OncoExperiment(project_name = "Example")
example <- load_assays(example, 
                       AssayTypes$RNA, 
                       data = "assets/rna_example.csv", 
                       metadata = "assets/rna_metadata_example.csv")
example
```

    ## An OncoExperiment object
    ## Project Name: Example 
    ## Version: 0.1.0 
    ## Experiments: RNA 
    ## 
    ## Experiment summary:
    ##  experiment    class    dim
    ##         RNA RNAAssay 5 x 49

#### Logical Indexing

You can use logical operations to subset data based on the rows or
columns using the `$` operator. Like in other libraries, you can use
`==`, `!=`, `<`, `>`, `<=`, and `>=` to filter your data. All
experiments and assays will be automatically harmonized.

In our example data, the metadata file contains a column called
`AgeCategory`. Let’s say we want to subset by all samples where patients
are adults and filter out unknown or pediatric samples. We can do this
by creating a logical vector and using it to subset the `OncoExperiment`
object:

``` r

is_metastatic <- example$PrimaryOrMetastasis == "Metastatic" # returns logical vector
metastatic <- example[, is_metastatic]
metastatic
```

    ## An OncoExperiment object
    ## Project Name: Example 
    ## Version: 0.1.0 
    ## Experiments: RNA 
    ## 
    ## Experiment summary:
    ##  experiment    class    dim
    ##         RNA RNAAssay 5 x 20

We removed ~30 samples and were left with only samples that met our
criterion. If we had additional assays, they would have been
simultaneously reduced.

#### Standard Indexing

If you require more fine-grained control, you can subset by specific
indices:

``` r
OncoExperiment[i = rownames, j = primary or colnames, k = assay]
```

For example, we can do the same example as above by using standard
indexing:

``` r

example[, example$PrimaryOrMetastasis == "Metastatic", "RNA"]
```

    ## An OncoExperiment object
    ## Project Name: Example 
    ## Version: 0.1.0 
    ## Experiments: RNA 
    ## 
    ## Experiment summary:
    ##  experiment    class    dim
    ##         RNA RNAAssay 5 x 20

#### Functions

`MultiAssayExperiment` has a number of useful functions for subsetting
or intersecting. The previous examples will likely be the most common
ones you’ll use. For more specific use cases, please refer to the
[MultiAssayExperiment
documentation](https://waldronlab.io/MultiAssayExperiment/reference/subsetBy.html).

## `*Assay`

A select number of data formats for different types of experiments are
currently supported and have been inspired by the oncogenic files hosted
at [DepMap](https://depmap.org/portal/). Currently, we support:

| Name | Purpose |
|:---|:---|
| [`RNAAssay`](https://tycgenomicslab.github.io/articles/RNAAssay.md) | Cell line vs. gene expression matrix for RNA-seq data. Useful for understanding gene expression patterns in different cell lines. |
| [`DrugResponseAssay`](https://tycgenomicslab.github.io/articles/DrugResponseAssay.md) | Drugs vs cell lines with log-fold change drug response values. Useful for [PRISM](https://depmap.org/repurposing/)-like data to understand therapeutic responses. |

### Structure

All assays have their own S4 class but ultimately inherit from
[`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html).
This means that you can use the same accessors and methods to interact
with them. For example, you can use
[`assays()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
to access the raw values,
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
to access any metadata associated with the rows, and
[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
to access any metadata associated with the columns. Here is a helpful
graphic:

> ![](assets/se.png)\> Credit: [Morgan M, Obenchain V, Hester J, Pagès H
> (2026). SummarizedExperiment: A container (S4 class) for matrix-like
> assays. doi:10.18129/B9.bioc.SummarizedExperiment. R package version
> 1.42.0,
> https://bioconductor.org/packages/SummarizedExperiment.](https://www.bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)

### Enumerator

This package provides a helpful enumerator for seeing and accessing all
available assay types. This just helps with code suggestions, prevents
typos, and makes it easy to see what’s available.

``` r

AssayTypes # see all options
```

    ## # A generic enum: 2 members
    ##  chr RNA          : RNA
    ##  chr DrugResponse : DrugResponse

``` r

AssayTypes$RNA # easy to access their character vector
```

    ## [1] "RNA"

``` r

AssayTypes$DoesNotExist # if you make a typo, you'll get an error
```

    ## Error in `error_undefined_member()`:
    ## ! Cannot subset an undefined or unknown member 'DoesNotExist'

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] dplyr_1.2.1                 MultiAssayExperiment_1.38.0
    ##  [3] SummarizedExperiment_1.42.0 Biobase_2.72.0             
    ##  [5] GenomicRanges_1.64.0        Seqinfo_1.2.0              
    ##  [7] IRanges_2.46.0              S4Vectors_0.50.2           
    ##  [9] BiocGenerics_0.58.1         generics_0.1.4             
    ## [11] MatrixGenerics_1.24.0       matrixStats_1.5.0          
    ## [13] OncoProfilingTools_0.1.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10          SparseArray_1.12.2   lattice_0.22-9      
    ##  [4] digest_0.6.39        magrittr_2.0.5       evaluate_1.0.5      
    ##  [7] grid_4.6.1           fastmap_1.2.0        jsonlite_2.0.0      
    ## [10] Matrix_1.7-5         enumr_0.0.1          textshaping_1.0.5   
    ## [13] jquerylib_0.1.4      abind_1.4-8          cli_3.6.6           
    ## [16] crayon_1.5.3         rlang_1.3.0          XVector_0.52.0      
    ## [19] withr_3.0.3          cachem_1.1.0         DelayedArray_0.38.2 
    ## [22] yaml_2.3.12          BiocBaseUtils_1.14.2 otel_0.2.0          
    ## [25] S4Arrays_1.12.0      tools_4.6.1          vctrs_0.7.3         
    ## [28] R6_2.6.1             lifecycle_1.0.5      fs_2.1.0            
    ## [31] ragg_1.5.2           pkgconfig_2.0.3      desc_1.4.3          
    ## [34] pkgdown_2.2.1        bslib_0.12.0         pillar_1.11.1       
    ## [37] data.table_1.18.6.1  glue_1.8.1           systemfonts_1.3.2   
    ## [40] xfun_0.60            tibble_3.3.1         tidyselect_1.2.1    
    ## [43] knitr_1.52           htmltools_0.5.9      rmarkdown_2.32      
    ## [46] compiler_4.6.1
