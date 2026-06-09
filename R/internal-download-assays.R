# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(httr2)
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

#' Entry URL's from the API are not valid download links.
#' Instead, we can use them to get a csv table listing all
#' files available for download to find the correct download URL
#'
#' You can provide an OncoExperiment object to store the data.
#' Caching is useful as it takes a few secodnds to download the data.
#' This is full optional.
#'
#' @param object An optional OncoExperiment object to store the downloaded assay data. Default is NULL.
#' @return A data frame of available assay names.
#' @keywords internal
fetch_all_assay_urls <- function(object = NULL) {
  # If object given, check type and see if it is already stored in @metadata
  if (!is.null(object)) {
    if (!inherits(object, "OncoExperiment")) {
      stop("`object` must be an OncoExperiment object.")
    }

    if ("depmap_assays" %in% names(object@metadata)) {
      message("Using cached assay URLs from object metadata.")
      return(object@metadata$depmap_assays)
    }
  }

  url <- "https://depmap.org/portal/api/download/files"

  response_text <- httr2::request(url) |>
    httr2::req_user_agent("OncoExperiment R package email authored by lapierreja@vcu.edu") |>
    httr2::req_perform() |>
    httr2::resp_body_string()

  files <- read.csv(text = response_text, stringsAsFactors = FALSE)

  # Extract release, filename, and url as a data frame
  files_df <- data.frame(
    release = files$release,
    filename = files$filename,
    url = files$url,
    stringsAsFactors = FALSE
  )

  # If object given, store the data frame in @metadata for future use
  if (!is.null(object)) {
    object@metadata$depmap_assays <- files_df
  }

  files_df
}

#' Find the download URL for a specific assay file based on the release and filename
#' @param release The DepMap release, e.g. "25Q3"
#' @param filename The name of the file to download, e.g. "depmap_25Q3_expression.csv"
#' @return The download URL for the specified file
#' @keywords internal
find_assay_url <- function(release, filename) {
  files <- fetch_all_assay_urls()

  # Find the row where "DepMap Public <release>" mathces for the release
  # and the filename matches for the filename
  release_string <- paste("DepMap Public", release, sep = " ")
  matching_row <- files[files$release == release_string & files$filename == filename, , drop = FALSE]

  if (nrow(matching_row) == 0) {
    stop(paste("No matching file found for release", release, "and filename", filename))
  }

  # Return the URL as a single string, subsetting keeps it a character vector
  matching_row$url[1]
}

#' @keywords internal
download_assay <- function(release, url, filename) {
  if (!nzchar(url)) {
    stop("`url` cannot be empty.")
  }

  if (!nzchar(release)) {
    stop("`release` cannot be empty.")
  }

  if (!nzchar(filename)) {
    stop("`filename` cannot be empty.")
  }

  cache_dir <- file.path(".cache", release)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  cached_file_path <- file.path(cache_dir, filename)
  if (file.exists(cached_file_path)) {
    message(paste("File already exists in cache:", cached_file_path))
    return(cached_file_path)
  }

  message(paste("Downloading", filename))

  response <- httr2::request(url) |>
    httr2::req_user_agent("OncoExperiment R package email authored by lapierreja@vcu.edu") |>
    httr2::req_progress() |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_perform()

  body_raw <- httr2::resp_body_raw(response)
  content_type <- tolower(httr2::resp_header(response, "content-type"))
  if (grepl("text/html|application/xhtml\\+xml", content_type)) {
    stop(paste("Download returned HTML instead of CSV for", filename))
  }

  writeBin(body_raw, con = cached_file_path)

  if (!file.exists(cached_file_path)) {
    stop(paste("Failed to download file:", cached_file_path))
  }

  message(paste("File downloaded and cached at:", cached_file_path))

  cached_file_path
}
