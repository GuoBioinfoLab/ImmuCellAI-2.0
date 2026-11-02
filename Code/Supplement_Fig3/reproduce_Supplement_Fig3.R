# Frozen-data analysis and R plotting; no deconvolution, model selection or simulation.
local({
source_files <- Filter(function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x),
                       lapply(sys.frames(), function(x) x$ofile))
if (length(source_files)) {
  entry <- tail(source_files, 1L)[[1L]]
  args <- getOption("Supplement_Fig3.args", "--plots")
} else {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Use source() or Rscript.")
  entry <- sub("^--file=", "", file_arg[1L])
  args <- commandArgs(TRUE)
}
code_dir <- dirname(normalizePath(entry, winslash = "/", mustWork = TRUE))
main <- function() {
  option <- function(flag, default) {
    i <- match(flag, args)
    if (is.na(i)) return(default)
    if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value for ", flag)
    args[i + 1L]
  }
  data_dir <- normalizePath(option("--data-dir", getOption("Supplement_Fig3.data_dir",
                          file.path(code_dir, "../../Data/Supplement_Fig3"))), winslash = "/", mustWork = TRUE)
  output_dir <- option("--output-dir", getOption("Supplement_Fig3.output_dir", file.path(code_dir, "output_plot")))
  panels <- option("--panel", getOption("Supplement_Fig3.panels", c("A", "B")))
  if (identical(panels, "all")) panels <- c("A", "B")
  output_format <- option("--format", getOption("Supplement_Fig3.format", "both"))
  dpi <- as.integer(option("--dpi", getOption("Supplement_Fig3.dpi", 300L)))
  stopifnot(is.character(args), length(panels) > 0L, !anyDuplicated(panels),
            all(panels %in% c("A", "B")), length(output_format) == 1L,
            output_format %in% c("pdf", "png", "both"), length(dpi) == 1L, !is.na(dpi), dpi >= 72)
  had_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", globalenv())
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, globalenv())
    else if (exists(".Random.seed", globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  source(file.path(code_dir, "shared/statistics.R"), local = TRUE, encoding = "UTF-8")
  context <- prepare_analysis(data_dir)
  analysis <- list()
  analysis$A <- list(points = context$scatter, metrics = context$metrics)
  analysis$B <- context$summary
  cat("PASS: Supplement_Fig3 frozen-data analysis and numeric checks.\n")
  files <- character()
  if (("--plots" %in% args || "--panel" %in% args) && !"--verify-only" %in% args) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
    if (tolower(output_dir) == tolower(data_dir) ||
        startsWith(paste0(tolower(output_dir), "/"), paste0(tolower(data_dir), "/")))
      stop("Output directory cannot be inside the Data directory.")
    source(file.path(code_dir, "plot_panels.R"), local = TRUE, encoding = "UTF-8")
    cat("PDF/PNG output:", output_dir, "\n")
  }
  invisible(list(analysis = analysis, verification = context$verification, files = files))
}
main()
})
