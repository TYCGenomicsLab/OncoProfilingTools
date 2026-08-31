# ExpressionAssay S4 class definition

\`ExpressionAssay\` is a \`SummarizedExperiment\` subclass for storing
gene or transcript expression data.

## Value

An object of class \`ExpressionAssay\`.

## Details

The primary assay is named \`"expression"\` and should contain a numeric
matrix-like object with features as rows and samples/models as columns.

This class follows the Bioconductor \`SummarizedExperiment\` convention:

\- rows represent features, such as genes or transcripts - columns
represent samples, models, or cell lines - \`rowData()\` stores
feature-level metadata - \`colData()\` stores sample/model-level
metadata - \`metadata()\` stores assay-level metadata such as units and
provenance
