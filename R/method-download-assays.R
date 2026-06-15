# Author: Jason LaPierre
# Last update: June 4, 2026
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)
library(dplyr)

#' Interacts with the DepMap API to download the specified supported assays.
#' It will download the assay files for the latest release of that type in store them in the `.cache` directory.
#'
#' If the files exist already, it will skip their download.
#' @name download_assays
#' @param object An OncoExperiment object to store the downloaded assay data.
#' @param assay_type The type of assay to download, e.g. "expression". See supported assay types with `list_supported_assays()`.
#' @include internal-load-assays.R
#' @export
setMethod("download_assays", "OncoExperiment", function(
  object,
  assay_type,
  ...
) {
  # if assay_type not in supported_assays, throw an error
  if (!assay_type %in% supported_assays$assay_type) {
    stop(paste("Unsupported assay type:", assay_type, "\nSupported assay types are:", paste(unique(supported_assays$assay_type), collapse = ", ")))
  }

  available_assays <- fetch_all_assay_urls(object)

  supported_types <- supported_assays |>
    dplyr::filter(is.na(assay_type) | assay_type == .env$assay_type)

  assay_urls <- available_assays |>
    dplyr::inner_join(
      supported_types,
      by = c("release" = "release_name_pattern"),
      relationship = "many-to-many"
    ) |>
    dplyr::rowwise() |>
    dplyr::filter(grepl(filename_pattern, filename)) |>
    dplyr::ungroup() |>
    dplyr::select(release, url, filename, metadata_assay) |>
    dplyr::distinct()

  if (nrow(assay_urls) == 0) {
    stop(paste("No supported DepMap files found for assay type", assay_type))
  }

  ## -----------------------------------------------
  ## Download each assay file and store in .cache directory
  ## TODO: implement custom directory option
  ## -----------------------------------------------
  for (selection in seq_len(nrow(assay_urls))) {
    download_assay(
      assay_urls$release[selection],
      assay_urls$url[selection],
      assay_urls$filename[selection]
    )
  }

  ## -----------------------------------------------
  ## If there are associated metadata files, let's
  ## recur and download those too
  ## -----------------------------------------------
  metadata_assays <- assay_urls$metadata_assay |> stats::na.omit() |> unique()
  for (metadata_assay in metadata_assays) {
    download_assays(object, assay_type = metadata_assay)
  }
})
