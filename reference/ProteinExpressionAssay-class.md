# ProteinExpressionAssay S4 class definition

\`ProteinExpressionAssay\` is a \`SummarizedExperiment\` subclass for
storing protein-coding gene expression data.

## Value

An object of class \`ProteinExpressionAssay\`.

## Details

The primary assay is named \`"protein_expression"\` and should contain a
numeric matrix-like object with features as rows and samples/models as
columns.

This class follows the Bioconductor \`SummarizedExperiment\` convention:

\- rows represent features, such as protein-coding genes - columns
represent samples, models, or cell lines - \`rowData()\` stores
feature-level metadata - \`colData()\` stores sample/model-level
metadata - \`metadata()\` stores assay-level metadata such as units and
provenance
