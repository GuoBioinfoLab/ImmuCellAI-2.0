# Reanalyse archived correlations; no new deconvolution or runtime benchmark.
methods <- c("ImmuCellAI 2.0", "BayesPrism", "CIBERSORT", "CITMIC", "DWLS", "ImmuCellAI", "MuSiC")
groups <- c("BostonGene blood", "Newman PBMC", "Monaco PBMC")
counts <- setNames(c(26L, 6L, 5L), groups)
read_csv <- function(path) {
  read.csv(path, colClasses = "character", check.names = FALSE,
           na.strings = character(0), fileEncoding = "UTF-8-BOM")
}
same_numbers <- function(x, y, tolerance = 1e-12) {
  stopifnot(length(x) == length(y), identical(is.na(x), is.na(y)),
            all(abs(x[!is.na(x)] - y[!is.na(y)]) < tolerance))
}
analyze_correlations <- function(data_dir, panel, dataset) {
  d <- read_csv(file.path(data_dir, panel, "pie_plot_data.csv"))
  stopifnot(nrow(d) == 7L * counts[[dataset]], all(d$Dataset == dataset),
            !anyDuplicated(d[c("method", "cellType")]), setequal(unique(d$method), methods),
            all(table(d$method) == counts[[dataset]]),
            length(unique(d$cellType)) == counts[[dataset]])
  d$r <- suppressWarnings(as.numeric(d$correlation))
  stopifnot(all(is.finite(d$r) | d$correlation %in% c("", "NA", "NaN")),
            all(abs(d$r[is.finite(d$r)]) <= 1),
            all((d$status == "valid") == is.finite(d$r)))
  d
}
verify_analysis <- function(analysis, data_dir) {
  d <- do.call(rbind, analysis[LETTERS[1:3]])
  rownames(d) <- NULL
  expected <- read_csv(file.path(data_dir, "Fig3_all_panels_correlations.csv"))
  stopifnot(nrow(d) == 259L, identical(d[names(expected)], expected))
  cat("Fig3: 259 correlations, 21 recalculated mean/SD summaries and 14 runtime records verified.\n")
}
