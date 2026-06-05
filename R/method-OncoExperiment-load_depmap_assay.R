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
#' @include internal-OncoExperiment-load_assays.R
#' @return An updated OncoExperiment object with the downloaded assay data stored in the @assays slot.
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
      by = c("release" = "release_name_pattern")
    ) |>
    dplyr::filter(grepl(filename_pattern, filename)) |>
    dplyr::select(release, url, filename) |>
    dplyr::distinct()

  if (nrow(assay_urls) == 0) {
    stop(paste("No supported DepMap files found for assay type", assay_type))
  }

  # Example url: /portal/data_page/?release=DepMap+Public+&file=FILENAME
  # Substitute the 26Q1 in each URL with the release parameter provided by the user.
  # It may not neccessarily be 26Q1, but it is the last item in the release parameter
  for (selection in seq_len(nrow(assay_urls))) {
    download_assay(
      assay_urls$release[selection],
      assay_urls$url[selection],
      assay_urls$filename[selection]
    )
  }

  ## For now, return the unchanged object so package sourcing succeeds.
  return(object)
})
