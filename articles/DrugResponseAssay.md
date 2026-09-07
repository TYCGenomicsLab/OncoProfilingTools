# Drug Response Assay (LFC)

## Format

For a log-fold change drug response assay (DRA), your data should be
formated with cell lines as columns and drugs as rows. The first column
is expected to be the drug identifiers. For example:

![](assets/dra_format.png)

### Metadata

Providing metadata is optional but will likely accompany your data.
Metadata should be provided in a separate tab-delimited file. The first
column should still be the same drug identifiers that can be mapped to
those listed in the data file. Here is an example of a metadata file:

| Arbitrary.ID | Drug.Name..Brand. | Drug.Name..Generic. | Drug.Company | Dosage |
|:---|:---|:---|:---|:---|
| DRG-001 | Tobramycin Sulfate | TOBRAMYCIN SULFATE |  | 9.1mg |
| DRG-002 | Aptivus | tipranavir | Boehringer Ingelheim Pharmaceuticals, Inc. | 6.0mg |
| DRG-003 | Ketoprofen | ketoprofen |  | 5.1mg |
| DRG-004 | Nortriptyline Hydrochloride | Nortriptyline Hydrochloride | REMEDYREPACK INC. | 9.4mg |
| DRG-005 | Diovan HCT | valsartan and hydrochlorothiazide | PD-Rx Pharmaceuticals, Inc. | 4.0mg |
| DRG-006 | napoleon PERDIS SHEER GENIUS LIQUID FOUNDATION BROAD SPECTRUM SPF 20 Look 5 | OCTINOXATE, TITANIUM DIOXIDE | Napoleon Perdis Cosmetics, Inc | 2.5mg |

## Loading

Loading a DRA is achieved similarly to other assays. For this example we
will load in metadata but this is optional. We’ve created `csv` files
with some 50x50 dummy data to act as our example:

``` r

library(OncoProfilingTools)
library(MultiAssayExperiment)
library(SummarizedExperiment)

exp <- OncoExperiment(project_name = "Example DRA")
exp <- load_assays(
  object = exp,
  assay_type = AssayTypes$DrugResponse,
  data = "assets/dra_example.csv",
  metadata = "assets/dra_metadata_example.csv"
)

exp
```

    ## An OncoExperiment object
    ## Project Name: Example DRA 
    ## Version: 0.1.0 
    ## Experiments: DrugResponse 
    ## 
    ## Experiment summary:
    ##    experiment             class     dim
    ##  DrugResponse DrugResponseAssay 50 x 50

## Navigating data

Your `OncoExperiment` links all your experimental assays using the
identifiers provided in your files. We can check that they loaded in
correctly by observing the `colData` (see
[MultiAssayExperiment](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/QuickStartMultiAssay.html)
for more information on this):

``` r

head(colData(exp))
```

    ## DataFrame with 6 rows and 1 column
    ##            model_id
    ##         <character>
    ## VQM-732     VQM-732
    ## 0Xe-933     0Xe-933
    ## dFL-450     dFL-450
    ## EcW-351     EcW-351
    ## nkn-474     nkn-474
    ## yBR-806     yBR-806

### Access raw values

If you need to look at the raw values of your data or interact with them
directly, you can access the `DrugResponseAssay` (instance of
`SummarizedExperiment`) object directly.

``` r

dra <- experiments(exp)$DrugResponse # access by whatever you named the assay
dra
```

    ## class: DrugResponseAssay 
    ## dim: 50 50 
    ## metadata(8): assay_name source ... n_missing_values
    ##   missing_value_percent
    ## assays(1): response
    ## rownames(50): DRG-001 DRG-002 ... DRG-049 DRG-050
    ## rowData names(5): Arbitrary.ID Drug.Name..Brand. Drug.Name..Generic.
    ##   Drug.Company Dosage
    ## colnames(50): VQM-732 0Xe-933 ... IuA-546 wde-847
    ## colData names(1): model_id

#### Get matrix of raw values

The raw values are stored in the `assays` slot of the
`SummarizedExperiment` object. Accessing them is relatively
straightforward. In this example, we’ll just print the first 3 rows of
the matrix as an example. You should see drug names as row names and
cell line names as column names.

``` r

dra_matrix <- assays(dra)$response # TODO: future function will make this easier (i.e. get_values(...)
knitr::kable(head(dra_matrix, n = 3L))
```

|  | VQM-732 | 0Xe-933 | dFL-450 | EcW-351 | nkn-474 | yBR-806 | 6rZ-921 | CVp-900 | FBr-666 | pwX-711 | 8Ak-004 | Wvk-040 | uuk-615 | wwg-515 | 1Qu-701 | JPn-336 | 5Dj-010 | lW5-428 | hGO-989 | eIt-266 | Bi0-209 | mMA-166 | iA6-055 | GRe-728 | ipU-508 | r6v-183 | OV3-829 | 6oP-633 | Vd7-714 | wyH-085 | ckl-460 | FdA-492 | gtA-835 | FT5-355 | fws-845 | aYK-906 | CyJ-883 | w15-391 | spU-324 | qK6-434 | Y4t-949 | jpJ-859 | IK1-886 | NuV-036 | oFd-671 | Aeb-277 | zSF-424 | BR9-197 | IuA-546 | wde-847 |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DRG-001 | 2.78854 | NA | NA | -5.53579 | NA | NA | 7.84359 | -8.26122 | -1.56156 | -9.40406 | NA | 0.10711 | -9.46928 | -6.02325 | 2.99769 | 0.89883 | -5.59119 | 1.78531 | 6.18861 | -9.87002 | 6.11639 | 3.96279 | -3.19499 | -6.89041 | 9.14426 | -3.26811 | NA | -8.06567 | 6.94989 | 2.07452 | 6.14257 | 4.59464 | NA | NA | NA | 1.04081 | 6.58809 | 2.37040 | NA | NA | NA | -9.08351 | -5.44203 | -4.21224 | NA | -5.34418 | NA | -4.44053 | NA | NA |
| DRG-002 | NA | NA | -4.66044 | 8.73309 | 2.96071 | 2.18262 | NA | 4.58254 | -6.73195 | -2.41089 | 9.79047 | 2.80000 | 1.13899 | 3.69229 | 6.85704 | 5.52000 | NA | -9.35800 | NA | -4.64518 | -5.78034 | 8.85819 | 7.52735 | -3.70644 | 3.10877 | -2.08736 | 8.29095 | NA | -4.70240 | NA | 1.22736 | -4.74517 | 1.69172 | 7.95646 | -2.01199 | -5.61358 | 9.95075 | 0.19053 | -8.18181 | NA | -7.80702 | 2.54892 | 5.84159 | -1.55680 | -8.72945 | -2.36761 | 9.92243 | NA | 9.42157 | 7.21559 |
| DRG-003 | -9.77038 | 4.41444 | 3.63421 | 0.73941 | -4.66350 | 2.81924 | NA | -1.30469 | -0.92553 | NA | 7.51706 | -4.73222 | 0.01172 | -6.42696 | 8.25256 | NA | -4.03110 | 2.77899 | 2.17940 | -6.94321 | NA | NA | 5.57253 | 0.60707 | -9.98856 | -3.51688 | -9.61047 | NA | 7.57444 | 6.63331 | -3.84972 | -8.84150 | 7.56019 | 8.93899 | -8.28693 | -0.28019 | -8.61575 | 5.21204 | 5.31669 | -7.43217 | -0.49435 | NA | -4.69887 | 7.44866 | NA | -5.76404 | 0.78592 | 4.59862 | -5.97698 | -3.76567 |

#### Look at metadata

If you loaded metadata alongside the assay, it will be stored in the
`rowData` slot of the `SummarizedExperiment` object. The `rownames` are
the drug identifiers and `rowData` is the associated metadata.

``` r

metadata <- rowData(dra)
knitr::kable(head(metadata, n = 3L))
```

|  | Arbitrary.ID | Drug.Name..Brand. | Drug.Name..Generic. | Drug.Company | Dosage |
|:---|:---|:---|:---|:---|:---|
| DRG-001 | DRG-001 | Tobramycin Sulfate | TOBRAMYCIN SULFATE |  | 9.1mg |
| DRG-002 | DRG-002 | Aptivus | tipranavir | Boehringer Ingelheim Pharmaceuticals, Inc. | 6.0mg |
| DRG-003 | DRG-003 | Ketoprofen | ketoprofen |  | 5.1mg |

> ⚠️ This is not to be confused with the `metadata()` function or slot.
> This returns a DataFrame of any metadata associated with the assay,
> not the samples themselves. This could contain anything from citation
> information, data about the file source, etc.
