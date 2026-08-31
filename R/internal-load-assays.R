# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)


#' Load a single assay using a registered loader
#'
#' @param path Character scalar. Path to the assay file.
#' @param filename Character scalar. File name being loaded.
#' @param assay_type Character scalar. Project-level assay type.
#' @param OncoAssay_class Character scalar. Expected returned assay class.
#' @param loader Character scalar. Registered loader function name.
#' @param ... Additional arguments passed to the registered loader.
#'
#' @return An assay object of the expected class.
#' @keywords internal
.load_assay <- function(
  assay_type,
  OncoAssay_class,
  loader,
  data_path,
  metadata_path = NULL,
  identifier = NULL,
  name = NULL,
  ...
) {
  if (!is.character(data_path) || length(data_path) != 1L || is.na(data_path) || !nzchar(data_path)) {
    stop("`data_path` must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(data_path)) {
    stop("Assay file does not exist: ", data_path, call. = FALSE)
  }

  if (!is.null(metadata_path)) {
    if (!is.character(metadata_path) || length(metadata_path) != 1L || is.na(metadata_path) || !nzchar(metadata_path)) {
      stop("`metadata_path` must be a single non-empty character string.", call. = FALSE)
    }
    if (!file.exists(metadata_path)) {
      stop("Metadata file does not exist: ", metadata_path, call. = FALSE)
    }
  }

  if (!is.character(assay_type) || length(assay_type) != 1L || is.na(assay_type) || !nzchar(assay_type)) {
    stop("`assay_type` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(OncoAssay_class) || length(OncoAssay_class) != 1L || is.na(OncoAssay_class) || !nzchar(OncoAssay_class)) {
    stop("`OncoAssay_class` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(loader) || length(loader) != 1L || is.na(loader) || !nzchar(loader)) {
    stop("`loader` must be a single non-empty character string.", call. = FALSE)
  }

  if (!exists(loader, mode = "function")) {
    stop(
      "Registered loader function does not exist: ",
      loader,
      call. = FALSE
    )
  }

  loader_fn <- get(loader, mode = "function")

  assay <- loader_fn(
    data_path = data_path,
    metadata_path = metadata_path,
    assay_name = assay_type,
    identifier = identifier,
    name = name,
    ...
  )

  if (!inherits(assay, OncoAssay_class)) {
    stop(
      "Loader `",
      loader,
      "` returned class `",
      class(assay)[1L],
      "`, but registry expected `",
      OncoAssay_class,
      "`.",
      call. = FALSE
    )
  }

  validObject(assay)

  assay
}
