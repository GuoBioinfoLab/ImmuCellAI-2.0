# Prepare local copies; never execute the archived analyses or modify inputs.
prepare_source_archive <- function(analysis_dir, prepared_dir, path_roots) {
  if (!is.character(path_roots) || is.null(names(path_roots)) ||
      anyDuplicated(names(path_roots)) || anyNA(path_roots) ||
      any(!nzchar(names(path_roots))) || any(!nzchar(path_roots))) {
    stop("path_roots must be a uniquely named, nonempty character vector.")
  }
  path_roots <- chartr("\\", "/", path_roots)
  if (any(!grepl("^(/|[A-Za-z]:/)", path_roots)) ||
      any(grepl("[\"'\r\n]", path_roots))) {
    stop("Use absolute paths without quotes or line breaks in path_roots.")
  }
  analysis_dir <- normalizePath(analysis_dir, winslash = "/", mustWork = TRUE)
  prepared_dir <- normalizePath(prepared_dir, winslash = "/", mustWork = FALSE)
  if (dir.exists(prepared_dir) || file.exists(prepared_dir)) {
    stop("Refusing to overwrite an existing prepared directory: ", prepared_dir)
  }
  files <- list.files(analysis_dir, pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE)
  files <- files[grepl("/source_archive/", files, fixed = TRUE)]
  if (!length(files)) stop("No archived R scripts found.")

  prepared <- lapply(files, function(f) {
    txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
    tokens <- unique(unlist(regmatches(txt, gregexpr("<LOCAL_[A-Z_]+>", txt))))
    required <- sub(">$", "", sub("^<", "", tokens))
    missing <- setdiff(required, names(path_roots))
    if (length(missing)) stop("Unconfigured roots in ", f, ": ", paste(missing, collapse = ", "))
    for (root in required) txt <- gsub(paste0("<", root, ">"), path_roots[[root]], txt, fixed = TRUE)
    parse(text = txt, keep.source = FALSE)
    relative <- substring(f, nchar(analysis_dir) + 2L)
    relative <- sub("/source_archive/", "/", relative, fixed = TRUE)
    list(path = file.path(prepared_dir, relative), text = txt)
  })
  # Validate every script before creating any destination files.
  for (item in prepared) {
    dir.create(dirname(item$path), recursive = TRUE, showWarnings = FALSE)
    writeLines(enc2utf8(item$text), item$path, useBytes = TRUE)
  }
  message("Prepared ", length(prepared), " scripts in ", prepared_dir,
          ". No analyses were run. See REPRODUCING_FIGURES.md for inputs and order.")
  invisible(vapply(prepared, `[[`, character(1), "path"))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: Rscript analysis/prepare_source_archive.R analysis/archive_paths.R")
  settings <- new.env(parent = baseenv())
  sys.source(args[1L], envir = settings)
  if (!exists("path_roots", settings, inherits = FALSE) ||
      !exists("prepared_dir", settings, inherits = FALSE)) {
    stop("Config must define path_roots and prepared_dir.")
  }
  prepare_source_archive("analysis", settings$prepared_dir, settings$path_roots)
}
