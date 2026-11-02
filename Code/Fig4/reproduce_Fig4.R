# Recalculate statistics from frozen Fig4 inputs, without changing the cohort,
# phenotype assignments, deconvolution outputs, or model predictions.
local({
# source() records the current filename in its call frame, not in --file.
source_files <- lapply(sys.frames(), function(frame) frame$ofile)
source_files <- Filter(function(path) is.character(path) && length(path) == 1L &&
                         !is.na(path) && nzchar(path), source_files)
if (length(source_files)) {
  script_file <- source_files[[length(source_files)]]
  args <- getOption("Fig4.args", "--plots")
} else {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(script_arg)) stop("Run this file with source() or Rscript.")
  script_file <- sub("^--file=", "", script_arg[1L])
  args <- commandArgs(trailingOnly = TRUE)
}
if (!is.character(args)) stop("Fig4.args must be a character vector.")
code_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
data_dir <- normalizePath(file.path(code_dir, "../../Data/Fig4"), winslash = "/")
out <- getOption("Fig4.output_dir", file.path(code_dir, "output_plot"))
if (!is.character(out) || length(out) != 1L || is.na(out) || !nzchar(out)) {
  stop("Fig4.output_dir must be one non-empty directory path.")
}
suppressPackageStartupMessages({library(survival); library(jsonlite); library(pROC)})
read_tab <- function(path) read.delim(file.path(data_dir, path), check.names = FALSE, stringsAsFactors = FALSE)
checks <- list()
check <- function(name, condition) {
  checks[[name]] <<- isTRUE(condition)
  if (!isTRUE(condition)) stop("Verification failed: ", name)
}
near <- function(x, y, tol = 1e-11) {
  x <- as.numeric(x); y <- as.numeric(y)
  length(x) == length(y) && identical(is.na(x),is.na(y)) &&
    all(abs(x[!is.na(x)]-y[!is.na(y)]) < tol)
}

# Preserve cross-panel checks while storing each panel in its own directory.
source(file.path(code_dir,"shared/statistics.R"), local=TRUE)
for (panel in c("A","B","C","D","F","G","H")) {
  source(file.path(code_dir,panel,"analyze.R"), local=TRUE)
}

cat("All", length(checks), "numerical checks passed. H AUC =", format(auc,digits=16), "\n")

if ("--plots" %in% args) source(file.path(code_dir,"plot_panels.R"), local=TRUE)
})
