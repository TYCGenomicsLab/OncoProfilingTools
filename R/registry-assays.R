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
  ~assay_name, ~OncoAssay_class, ~loader,

  ## ---------------------------------------------------
  ## PRISM response assay
  ## ---------------------------------------------------
  "DrugResponse", "DrugResponseAssay", "load_drug_response_assay",

  ## ---------------------------------------------------
  ## RNA expression assay
  ## ---------------------------------------------------
  "RNA", "RNAAssay", "load_RNA_assay"
)

#' List supported assay types
#'
#' Returns all supported assay types declared in the assay registry.
#'
#' @return A character vector of supported assay types.
#'
#' @export
get_supported_types <- function() {
  unique(stats::na.omit(supported_assays$assay_name))
}

#' List all supported assays
#'
#' Returns all supported assays in a tribble to be used as a lookup table.
#'
#' @return A tribble
#'
#' @export
get_supported_assays <- function() {
  supported_assays
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
    !is.na(supported_assays$assay_name) &
      !is.na(supported_assays$OncoAssay_class) &
      !is.na(supported_assays$loader), ,
    drop = FALSE
  ]

  if (!is.null(assay_type)) {
    registry <- registry[
      registry$assay_name %in% assay_type, ,
      drop = FALSE
    ]
  }

  registry
}

#' An enumeration of supported assay types
#'
#' Helpful for intelligent code completion and avoiding typos in assay types.
#'
#' @examples
#' AssayTypes$RNA
#' AssayTypes$PRISM
#'
#' @export
AssayTypes <- enumr::enum( # nolint
  RNA = "RNA",
  DrugResponse = "DrugResponse"
)
