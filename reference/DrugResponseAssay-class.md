# DrugResponseAssay S4 class definition

\`DrugResponseAssay\` is a \`SummarizedExperiment\` subclass for storing
drug response data.

## Value

An object of class \`DrugResponseAssay\`.

## Details

The primary assay is named \`"response"\` and should contain a numeric
matrix-like object with treatments/drugs as rows and samples/models as
columns.

This class follows the Bioconductor \`SummarizedExperiment\` convention:

\- rows represent treatments, drugs, or compounds - columns represent
samples, models, or cell lines - \`rowData()\` stores
treatment/compound-level metadata - \`colData()\` stores
sample/model-level metadata - \`metadata()\` stores assay-level metadata
such as units and provenance

Missing values are allowed because drug response matrices are often
partially observed.
