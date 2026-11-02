# Compatibility entry; analysis and plotting are organized by panel.
local({
  run <- function() {
  paths <- Filter(function(x) is.character(x) && length(x) == 1L &&
                    !is.na(x) && nzchar(x), lapply(sys.frames(), function(f) f$ofile))
  if (length(paths)) entry <- paths[[length(paths)]]
  else {
    entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(entry) != 1L) stop("Run this script with source() or Rscript.")
    previous <- getOption("Fig2.args")
    args <- commandArgs(TRUE)
    options(structure(list(c("--plots", args)), names = "Fig2.args"))
    on.exit(options(structure(list(previous), names = "Fig2.args")), add = TRUE)
  }
  path <- file.path(dirname(normalizePath(entry, winslash = "/", mustWork = TRUE)),
                    "reproduce_Fig2.R")
  invisible(source(path, local = TRUE, encoding = "UTF-8")$value)
  }
  run()
})
