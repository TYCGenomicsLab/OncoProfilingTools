library(tibble)

#' Mapping of supported assay types, their corresponding OncoAssay classes, and filename patterns
#' These are the only files downloaded from DepMap, inferred, and loadable in their assay-specific classes.
#' @keywords internal
supported_assays <- tribble(
  ~release_name_pattern, ~assay_type, ~OncoAssay_class, ~filename_pattern,
  "DepMap Public 26Q1", NA, NA, "README.txt", # support downloading README file
  "DepMap Public 26Q1", NA, NA, "Model.csv",
  "DepMap Public 26Q1", NA, NA, "Gene.csv",
  "DepMap Public 26Q1", "Expression", "ExpressionAssay", "OmicsExpressionTranscriptTPMLogp1HumanAllGenes.csv",
  "DepMap Public 26Q1", "Expression", "ExpressionAssay", "OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded.csv",
  "DepMap Public 26Q1", "Protein Expression", "ProteinExpressionAssay", "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  "DepMap Public 26Q1", "Protein Expression", "ProteinExpressionAssay", "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv",
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", NA, ".*_Cell_Line_Meta_Data\\.csv$",
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", NA, ".*_Treatment_Meta_Data\\.csv$",
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", NA, ".*_Extended_Primary_Compound_List\\.csv$",
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", NA, ".*Readme\\txt$",
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", "TreatmentAssay", ".*_Extended_Primary_Data_Matrix\\.csv$"
)