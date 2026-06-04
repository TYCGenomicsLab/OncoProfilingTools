# Author: Jason LaPierre
# Last update: June 4, 2026
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(dplyr)

#' downloadAssays method for OncoExperiment objects
#' Interacts with the DepMap API to download all assays for a specified release
#' @name load_depmap_assay
#' @export
setMethod("download_depmap_assays", "OncoExperiment", function(
  object,
  release = "26Q1"
) {
  requested_release <- paste0("DepMap Public ", release)

  download_assay <- function(release, filename, url) {
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

    response <- httr2::request(url) |>
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

    cached_file_path
  }

  # We'll fetch all current assays and URLs from the DepMap API and extract the URLS for requested data type
  assay_information <- fetch_all_assay_urls()

  # Data frame where data type matches the requested assay and download_entry_url has the URL
  assay_urls <- assay_information %>%
    filter(.data$release == requested_release) %>%
    select("release", "filename", "url")

  if (nrow(assay_urls) == 0) {
    stop(paste("No DepMap files found for release", release))
  }

  # Example url: /portal/data_page/?release=DepMap+Public+&file=FILENAME
  # Substitute the 26Q1 in each URL with the release parameter provided by the user.
  # It may not neccessarily be 26Q1, but it is the last item in the release parameter
  for (selection in seq_len(nrow(assay_urls))) {
    download_assay(
      assay_urls[["release"]][selection],
      assay_urls[["filename"]][selection],
      assay_urls[["url"]][selection]
    )
  }

  ## For now, return the unchanged object so package sourcing succeeds.
  return(object)
})
