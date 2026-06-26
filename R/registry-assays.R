# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(tibble)
library(enumr)

#' Mapping of supported assay types, their corresponding assay classes, filename patterns, and metadata dependencies.
#'
#' These are the only files downloaded from DepMap, inferred, and loadable in
#' their assay-specific classes.
#'
#' @keywords internal
supported_assays <- tibble::tribble(
  ~release_name_pattern,
  ~assay_type,
  ~OncoAssay_class,
  ~filename_pattern,
  ~loader,
  ~metadata_assay,

  ## ---------------------------------------------------
  ## DepMap Public 26Q1 metadata / auxiliary files
  ## ---------------------------------------------------

  "DepMap Public 26Q1",
  NA,
  NA,
  "README.txt",
  NA,
  NA,
  "DepMap Public 26Q1",
  "ExpressionMetadata",
  NA,
  "Model.csv",
  NA,
  NA,

  ## ---------------------------------------------------
  ## DepMap expression assays
  ## ---------------------------------------------------
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

  ## ---------------------------------------------------
  ## DepMap mutation assays
  ## ---------------------------------------------------

  "DepMap Public 26Q1",
  "Mutation",
  "MutationAssay",
  "OmicsSomaticMutations.csv",
  "load_mutation_assay",
  "ExpressionMetadata",

  ## ---------------------------------------------------
  ## PRISM metadata
  ## ---------------------------------------------------

  "PRISM Primary Repurposing DepMap Public 24Q2",
  "DrugResponseAssayMetadata",
  NA,
  ".*_Cell_Line_Meta_Data\\.csv$",
  NA,
  NA,
  "PRISM Primary Repurposing DepMap Public 24Q2",
  "DrugResponseAssayMetadata",
  NA,
  ".*_Treatment_Meta_Data\\.csv$",
  NA,
  NA,
  "PRISM Primary Repurposing DepMap Public 24Q2",
  "DrugResponseAssayMetadata",
  NA,
  ".*_Extended_Primary_Compound_List\\.csv$",
  NA,
  NA,
  "PRISM Primary Repurposing DepMap Public 24Q2",
  "DrugResponseAssayMetadata",
  NA,
  ".*README\\.txt$",
  NA,
  NA,

  ## ---------------------------------------------------
  ## PRISM response assay
  ## ---------------------------------------------------

  "PRISM Primary Repurposing DepMap Public 24Q2",
  "PRISM",
  "DrugResponseAssay",
  ".*_Extended_Primary_Data_Matrix\\.csv$",
  "load_drug_response_assay",
  "DrugResponseAssayMetadata"
)

#' List supported assay types
#'
#' Returns all supported assay types declared in the assay registry.
#'
#' @return A character vector of supported assay types.
#'
#' @export
get_supported_types <- function() {
  unique(stats::na.omit(supported_assays$assay_type))
}

#' Get loadable assay registry
#'
#' Returns a subset of the assay registry containing only rows that can be loaded
#' into assay objects.
#'
#' @param assay_type Optional character vector of assay types to keep.
#'
#' @return A tibble containing loadable assay registry rows.
#'
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

#' An enumeration of supported assay types
#'
#' Helpful for intelligent code completion and avoiding typos in assay type
#' strings.
#'
#' @examples
#' AssayTypes$Expression
#' AssayTypes$ProteinExpression
#' AssayTypes$PRISM
#'
#' @export
AssayTypes <- enumr::enum( # nolint
  ExpressionMetadata = "ExpressionMetadata",
  Expression = "Expression",
  ProteinExpression = "Protein Expression",
  DrugResponseAssayMetadata = "DrugResponseAssayMetadata",
  PRISM = "PRISM",
  Mutation = "Mutation",
  MutationMetadata = "MutationMetadata"
)
