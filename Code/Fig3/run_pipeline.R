# source() loads the pipeline API; it does not run inference or write files.
local({
  paths <- Filter(function(x) is.character(x) && length(x) == 1L &&
    !is.na(x) && nzchar(x), lapply(sys.frames(), function(f) f$ofile))
  if (length(paths)) entry <- paths[[length(paths)]] else {
    entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(entry) != 1L) stop("Run with source() or Rscript.")
  }
  code_dir <- dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))
  project <- normalizePath(file.path(code_dir, "../.."), winslash = "/", mustWork = TRUE)
  figure_name <- "Fig3"
  source(file.path(code_dir, "shared/pipeline.R"), local = TRUE, encoding = "UTF-8")
  invisible(list(inputs = pipeline_inputs, predict = pipeline_predict,
    statistics = pipeline_statistics, plot = pipeline_plot))
})
