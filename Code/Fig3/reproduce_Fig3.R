# Run with source() or Rscript; all analysis results stay in memory.
local({
  paths <- Filter(function(x) is.character(x) && length(x) == 1L &&
                    !is.na(x) && nzchar(x), lapply(sys.frames(), function(f) f$ofile))
  if (length(paths)) {
    entry <- paths[[length(paths)]]
    args <- getOption("Fig3.args", "--plots")
  } else {
    entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(entry) != 1L) stop("Run this script with source() or Rscript.")
    args <- commandArgs(TRUE)
  }
  code_dir <- dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))
  main <- function() {
    if (!is.character(args)) stop("Fig3.args must be a character vector.")
    settings <- list(
      data_dir = getOption("Fig3.data_dir", file.path(code_dir, "../../Data/Fig3")),
      output_dir = getOption("Fig3.output_dir", file.path(code_dir, "output_plot")),
      panels = getOption("Fig3.panels", LETTERS[1:5]),
      format = getOption("Fig3.format", "both"), dpi = getOption("Fig3.dpi", 300L))
    draw <- any(args %in% c("--plots", "--png", "--panel")) || any(startsWith(args, "--panel="))
    i <- 1L
    while (i <= length(args)) {
      flag <- sub("=.*$", "", args[i])
      keys <- c("--data-dir" = "data_dir", "--output-dir" = "output_dir",
                "--panel" = "panels", "--format" = "format", "--dpi" = "dpi")
      if (flag %in% names(keys)) {
        if (grepl("=", args[i], fixed = TRUE)) value <- sub("^[^=]*=", "", args[i])
        else {
          i <- i + 1L
          if (i > length(args) || startsWith(args[i], "--")) stop("Missing value for ", flag)
          value <- args[i]
        }
        settings[[keys[[flag]]]] <- value
      } else if (!flag %in% c("--plots", "--png", "--check-only", "--verify-only")) {
        stop("Unknown argument: ", args[i])
      }
      i <- i + 1L
    }
    if (any(args %in% c("--check-only", "--verify-only"))) draw <- FALSE
    data_dir <- normalizePath(settings$data_dir, winslash = "/", mustWork = TRUE)
    panels <- settings$panels
    output_format <- settings$format
    dpi <- as.integer(settings$dpi)
    stopifnot(is.character(panels), length(panels) > 0L, !anyDuplicated(panels),
              all(panels %in% LETTERS[1:5]), length(output_format) == 1L,
              output_format %in% c("pdf", "png", "both"), length(dpi) == 1L,
              is.finite(dpi), dpi >= 72L)
    source(file.path(code_dir, "shared/statistics.R"), local = TRUE, encoding = "UTF-8")
    analysis <- list()
    analysis$A <- analyze_correlations(data_dir, "A", "BostonGene blood")
    analysis$B <- analyze_correlations(data_dir, "B", "Newman PBMC")
    analysis$C <- analyze_correlations(data_dir, "C", "Monaco PBMC")
    for (panel in c("D", "E")) {
      source(file.path(code_dir, panel, "analyze.R"), local = TRUE, encoding = "UTF-8")
    }
    verify_analysis(analysis, data_dir)
    plot_files <- character()
    if (draw) {
      if (!capabilities("cairo")) stop("Cairo support is required for PDF/PNG export.")
      output_dir <- settings$output_dir
      stopifnot(is.character(output_dir), length(output_dir) == 1L,
                !is.na(output_dir), nzchar(output_dir))
      output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
      if (tolower(output_dir) == tolower(data_dir) ||
          startsWith(paste0(tolower(output_dir), "/"), paste0(tolower(data_dir), "/"))) {
        stop("Output cannot be inside the figure's Data directory.")
      }
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
      source(file.path(code_dir, "plot_panels.R"), local = TRUE, encoding = "UTF-8")
    }
    invisible(list(analysis = analysis, files = plot_files))
  }
  main()
})
