# Reproduce Fig6 from archived cell fractions; no upstream deconvolution.
local({
source_files <- lapply(sys.frames(), function(frame) frame$ofile)
source_files <- Filter(function(path) is.character(path) && length(path) == 1L &&
                         !is.na(path) && nzchar(path), source_files)
if (length(source_files)) {
  entry <- source_files[[length(source_files)]]
  args <- getOption("Fig6.args", "--plots")
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Run this file with source() or Rscript.")
  entry <- sub("^--file=", "", file_arg[1L])
  args <- commandArgs(trailingOnly = TRUE)
}
code_dir <- dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))

main <- function() {
  if (!is.character(args)) stop("Fig6.args must be a character vector.")
  option <- function(flag, default) {
    index <- match(flag, args)
    if (is.na(index)) return(default)
    if (index == length(args) || startsWith(args[index + 1L], "--")) stop("Value missing for ", flag)
    args[index + 1L]
  }
  data_dir <- option("--data-dir", getOption("Fig6.data_dir", file.path(code_dir, "../../Data/Fig6")))
  output_dir <- option("--output-dir", getOption("Fig6.output_dir", file.path(code_dir, "output_plot")))
  valid_path <- function(path) is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)
  if (!valid_path(data_dir) || !valid_path(output_dir)) stop("Input/output directories must be non-empty paths.")
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  dpi <- as.integer(option("--dpi", getOption("Fig6.dpi", 300L)))
  output_format <- option("--format", getOption("Fig6.format", "both"))
  panels <- option("--panel", getOption("Fig6.panels", LETTERS[1:5]))
  stopifnot(length(dpi) == 1L, !is.na(dpi), dpi >= 72,
            length(output_format) == 1L, output_format %in% c("both", "pdf", "png"),
            is.character(panels), length(panels) > 0L, !anyDuplicated(panels),
            all(panels %in% LETTERS[1:5]))

  source(file.path(code_dir, "shared/statistics.R"), local = TRUE, encoding = "UTF-8")
  inputs <- verify_inputs(data_dir)
  analysis <- list()
  for (panel in LETTERS[1:4]) {
    analysis[[panel]] <- analyze_box_panel(panel, inputs, data_dir)
  }
  source(file.path(code_dir, "E/analyze.R"), local = TRUE, encoding = "UTF-8")
  cat("PASS:", comparisons, "numeric comparisons; TB 92 runs, PRJNA683803 211 runs; 53 cell types.\n")
  cat("Maximum absolute numeric difference:", format(largest_difference, digits = 16), "\n")
  plot_files <- character()
  if ("--plots" %in% args || "--panel" %in% args) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) saved_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
      if (had_seed) assign(".Random.seed", saved_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
    if (tolower(output_dir) == tolower(data_dir) ||
        startsWith(paste0(tolower(output_dir), "/"), paste0(tolower(data_dir), "/"))) {
      stop("Output directory cannot be inside Data/Fig6.")
    }
    source(file.path(code_dir, "plot_panels.R"), local = TRUE, encoding = "UTF-8")
  }
  invisible(list(numeric_comparisons = comparisons, max_abs_numeric_difference = largest_difference,
                 panels = lapply(analysis, "[[", "report"), files = plot_files,
                 independent_statistics = inputs$recomputed_stats))
}
main()
})
