# Reproduce Fig5 from frozen cell fractions; no upstream deconvolution or selection.
local({
source_files <- lapply(sys.frames(), function(frame) frame$ofile)
source_files <- Filter(function(path) is.character(path) && length(path) == 1L &&
                         !is.na(path) && nzchar(path), source_files)
if (length(source_files)) {
  entry <- source_files[[length(source_files)]]
  args <- getOption("Fig5.args", "--plots")
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Run this file with source() or Rscript.")
  entry <- sub("^--file=", "", file_arg[1L])
  args <- commandArgs(trailingOnly = TRUE)
}
code_dir <- dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))

main <- function() {
  if (!is.character(args)) stop("Fig5.args must be a character vector.")
  option <- function(flag, default) {
    idx <- match(flag, args)
    if (is.na(idx)) return(default)
    if (idx == length(args) || startsWith(args[idx + 1L], "--")) stop("Value missing for ", flag)
    args[idx + 1L]
  }
  data_dir <- option("--data-dir", getOption("Fig5.data_dir", file.path(code_dir, "../../Data/Fig5")))
  output_dir <- option("--output-dir", getOption("Fig5.output_dir", file.path(code_dir, "output_plot")))
  valid_path <- function(path) is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)
  if (!valid_path(data_dir) || !valid_path(output_dir)) stop("Input/output directories must be non-empty paths.")
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  dpi <- as.integer(option("--dpi", getOption("Fig5.dpi", 300L)))
  output_format <- option("--format", getOption("Fig5.format", "both"))
  panels <- option("--panel", getOption("Fig5.panels", LETTERS[1:5]))
  stopifnot(length(dpi) == 1L, !is.na(dpi), dpi >= 72,
            length(output_format) == 1L, output_format %in% c("both", "pdf", "png"),
            is.character(panels), length(panels) > 0L, !anyDuplicated(panels),
            all(panels %in% LETTERS[1:5]))

  source(file.path(code_dir, "shared/statistics.R"), local = TRUE, encoding = "UTF-8")
  inputs <- verify_inputs(data_dir)
  fraction <- inputs$fraction
  meta <- inputs$meta
  flags <- inputs$flags
  analysis <- list()
  panel_counts <- list()
  for (panel in LETTERS[1:3]) {
    analysis[[panel]] <- analyze_box_panel(panel, flags, data_dir)
  }
  for (panel in c("D", "E")) {
    source(file.path(code_dir, panel, "analyze.R"), local = TRUE, encoding = "UTF-8")
  }
  cat("PASS:", length(numeric_errors), "numeric comparisons; 327 samples, 53 cell types.\n")
  cat("Maximum absolute numeric difference:", format(max(unlist(numeric_errors)), digits = 16), "\n")
  plot_files <- character()
  if ("--plots" %in% args || "--panel" %in% args) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) saved_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
      if (had_seed) assign(".Random.seed", saved_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    source(file.path(code_dir, "plot_panels.R"), local = TRUE, encoding = "UTF-8")
  }
  invisible(list(numeric_errors = numeric_errors, panel_counts = panel_counts, files = plot_files))
}
main()
})
