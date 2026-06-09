# Author: Jason LaPierre
# Virginia Commonwealth University, Katarzyna Tyc Lab


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

  files_to_load <- candidate_files[
    vapply(
      candidate_files,
      function(file) {
        .matches_loadable_assay_file(
          file = file,
          registry = registry
        )
      },
      logical(1)
    )
  ]

  if (length(files_to_load) == 0L) {
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

  unique(normalizePath(files_to_load, winslash = "/", mustWork = TRUE))
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
