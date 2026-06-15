# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(tibble)
library(enumr)

#' Mapping of supported assay types, their corresponding OncoAssay classes, and filename patterns
#' These are the only files downloaded from DepMap, inferred, and loadable in their assay-specific classes.
#' @keywords internal
supported_assays <- tribble(
  ~release_name_pattern, # Name matching release name in DepMap API portal
  ~assay_type, # Custom name for the type of assay, e.g. "Expression", "PRISM", etc.
  ~OncoAssay_class, # The S4 class that the assay will be loaded into, e.g. "ExpressionAssay". NA if not loadable.
  ~filename_pattern, # Pattern to match the filename of the assay file
  ~loader, # Function to load the assay file
  ~metadata_assay, # The type of metadata assay associated with this assay. Will be co- downloaded/loaded

  "DepMap Public 26Q1", NA, NA, "README.txt", NA, NA, # support downloading README file

  "DepMap Public 26Q1", "ExpressionMetadata", NA, "Model.csv", NA, NA,
  "DepMap Public 26Q1", "ExpressionMetadata", NA, "Gene.csv", NA, NA,
  "DepMap Public 26Q1",
  "Expression",
  "ExpressionAssay",
  "OmicsExpressionTranscriptTPMLogp1HumanAllGenes.csv",
  "load_expression_assay",
  "ExpressionMetadata",
  "DepMap Public 26Q1",
  "Protein Expression",
  "ProteinExpressionAssay",
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  "load_protein_expression_assay",
  "ExpressionMetadata",
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", NA, ".*_Cell_Line_Meta_Data\\.csv$", NA, NA,
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", NA, ".*_Treatment_Meta_Data\\.csv$", NA, NA,
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", NA, ".*_Extended_Primary_Compound_List\\.csv$", NA, NA,
  "PRISM Primary Repurposing DepMap Public 24Q2", "PRISM", NA, ".*Readme\\txt$", NA, NA,
  "PRISM Primary Repurposing DepMap Public 24Q2",
  "PRISM", "TreatmentAssay",
  ".*_Extended_Primary_Data_Matrix\\.csv$",
  NA,
  NA
)

#' A helper function to list all supported assay types that OncoExperiment can load.
#' @name get_supported_types
#' @export
get_supported_types <- function() {
  unique(stats::na.omit(supported_assays$assay_type))
}

# Helper function to get the registry of loadable assays, optionally filtered by assay type.
#' @keywords internal
.get_loadable_assay_registry <- function(assay_type = NULL) {
  registry <- supported_assays[
    !is.na(supported_assays$assay_type) &
      !is.na(supported_assays$OncoAssay_class) &
      !is.na(supported_assays$loader), ,
    drop = FALSE
  ]

  if (!is.null(assay_type)) {
    registry <- registry[
      registry$assay_type %in% assay_type, ,
      drop = FALSE
    ]
  }

  registry
}

#' An enumeration of suppported assay types.
#' Helpful for intelligent code completion and to avoid typos in assay type strings.
#'
#' Curious what types are available? Just print `AssayTypes`.
#' @examples
#' AssayTypes$Expression
#' AssayTypes$ProteinExpression
#' @export
AssayTypes <- enumr::enum( # nolint
  ExpressionMetadata = "ExpressionMetadata",
  Expression = "Expression",
  ProteinExpression = "Protein Expression",
  PRISM = "PRISM"
)
