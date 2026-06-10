# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab

library(methods)

#' Returns a tibble subset of the loadable assay registry based on the specified assay type.
#' @param path A character string specifying the directory or specific file path where the assay data is located. The function will look for files in this location that match the expected naming conventions for the specified assay type.
#' @param assay_type A character string specifying the type of assay being loaded (e.g., "Expression", "PRISM", etc.). This will determine the specific OncoAssay used and correct loading pattern. If NULL, all assay types in the registry will be considered.
#' @return A tibble containing the loadable assay registry entries that match the specified assay type.
#' @keywords internal
.resolve_files_to_load <- function(path, assay_type = NULL) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop(
      "`path` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  if (!file.exists(path)) {
    stop(
      sprintf(
        "Specified path '%s' does not exist. Run `download_assays()` first or provide `path` manually.",
        path
      ),
      call. = FALSE
    )
  }

  registry <- .get_loadable_assay_registry(assay_type = assay_type)

  if (nrow(registry) == 0L) {
    if (is.null(assay_type)) {
      stop(
        "No loadable assay types are currently registered.",
        call. = FALSE
      )
    }

    stop(
      "No loadable assay types matched: ",
      paste(assay_type, collapse = ", "),
      "\nCurrently loadable assay types are: ",
      paste(get_supported_types(), collapse = ", "),
      call. = FALSE
    )
  }

  candidate_files <- if (dir.exists(path)) {
    list.files(
      path = path,
      recursive = TRUE,
      full.names = TRUE
    )
  } else {
    path
  }

  candidate_files <- candidate_files[file.exists(candidate_files)]
  candidate_files <- candidate_files[!dir.exists(candidate_files)]

  if (length(candidate_files) == 0L) {
    stop(
      sprintf("No files were found under path '%s'.", path),
      call. = FALSE
    )
  }

  resolved_files <- lapply(candidate_files, function(file) {
    .match_loadable_assay_file(
      file = file,
      registry = registry
    )
  })

  resolved_files <- Filter(Negate(is.null), resolved_files)

  if (length(resolved_files) == 0L) {
    msg <- sprintf(
      "No loadable assay files were found under path '%s'.",
      path
    )

    if (!is.null(assay_type)) {
      msg <- paste0(
        msg,
        "\nRequested assay type(s): ",
        paste(assay_type, collapse = ", ")
      )
    }

    msg <- paste0(
      msg,
      "\nCurrently loadable assay types are: ",
      paste(get_supported_types(), collapse = ", ")
    )

    stop(msg, call. = FALSE)
  }

  resolved_files <- do.call(rbind, resolved_files)

  tibble::as_tibble(resolved_files)
}

#' @keywords internal
.match_loadable_assay_file <- function(file, registry) {
  filename <- basename(file)

  matched_registry <- registry[
    vapply(
      registry$filename_pattern,
      function(pattern) {
        grepl(pattern, filename)
      },
      logical(1)
    ),
    ,
    drop = FALSE
  ]

  if (nrow(matched_registry) == 0L) {
    return(NULL)
  }

  if (nrow(matched_registry) > 1L) {
    stop(
      "File matched multiple loadable assay registry entries: ",
      filename,
      "\nMatched assay types: ",
      paste(unique(matched_registry$assay_type), collapse = ", "),
      "\nMatched loaders: ",
      paste(unique(matched_registry$loader), collapse = ", "),
      call. = FALSE
    )
  }

  tibble::tibble(
    path = normalizePath(file, winslash = "/", mustWork = TRUE),
    filename = filename,
    release_name_pattern = matched_registry$release_name_pattern,
    assay_type = matched_registry$assay_type,
    OncoAssay_class = matched_registry$OncoAssay_class,
    filename_pattern = matched_registry$filename_pattern,
    loader = matched_registry$loader
  )
}

#' @keywords internal
.matches_loadable_assay_file <- function(file, registry) {
  filename <- basename(file)

  any(vapply(
    registry$filename_pattern,
    function(pattern) {
      grepl(pattern, filename)
    },
    logical(1)
  ))
}

#' @keywords internal
.load_assay <- function(
  path,
  filename,
  assay_type,
  OncoAssay_class,
  loader,
  ...
) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop("Assay file does not exist: ", path, call. = FALSE)
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

  loader_fn <- get(loader, mode = "function") # gets the correct function for loading that filenames assay

  assay <- loader_fn(path = path, ...)

  if (!inherits(assay, "OncoAssay")) {
    stop(
      "Loader `", loader, "` did not return an object inheriting from `OncoAssay`.",
      call. = FALSE
    )
  }

  if (!inherits(assay, OncoAssay_class)) {
    stop("Loader `", loader, "` returned class `", class(assay)[1L],
         "`, but registry expected `", OncoAssay_class, "`.",
         call. = FALSE
    )
  }

  validObject(assay)

  assay
}
