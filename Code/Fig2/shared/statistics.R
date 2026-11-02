# Frozen correlation tables are inputs, not sample-level prediction/truth pairs.
methods <- c("ImmuCellAI2", "BayesPrism", "CIBERSORT", "CITMIC",
             "DWLS", "ImmuCellAI", "MuSiC_basic")
datasets <- c("GSE164522", "GSE146771", "GSE176078", "GSE162025")
read_table <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             fileEncoding = "UTF-8-BOM", na.strings = c("", "NA", "NaN", "nan"))
}
same_numbers <- function(x, y, tolerance = 1e-12) {
  stopifnot(length(x) == length(y), identical(is.na(x), is.na(y)),
            all(abs(x[!is.na(x)] - y[!is.na(y)]) < tolerance))
}
check_rows <- function(actual, reference, keys, fields) {
  key <- function(d) do.call(paste, c(d[keys], sep = "\r"))
  a <- key(actual); b <- key(reference)
  stopifnot(!anyDuplicated(a), !anyDuplicated(b), setequal(a, b))
  reference <- reference[match(a, b), , drop = FALSE]
  for (field in fields) {
    if (is.numeric(actual[[field]])) same_numbers(actual[[field]], reference[[field]])
    else stopifnot(identical(as.character(actual[[field]]), as.character(reference[[field]])))
  }
}
analyze_correlations <- function(data_dir, panel, dataset, n_targets) {
  d <- read_table(file.path(data_dir, panel, "pie_plot_data.tsv"))
  stopifnot(nrow(d) == 7L * n_targets, all(d$dataset == dataset),
            all(d$n == 40L), length(unique(d$target)) == n_targets,
            !anyDuplicated(d[c("method", "target")]), setequal(unique(d$method), methods),
            all(d$status %in% c("valid", "unavailable", "constant_prediction", "undefined_correlation")),
            all((d$status == "valid") == is.finite(d$Pearson)),
            all(abs(d$Pearson[is.finite(d$Pearson)]) <= 1))
  source <- read_table(file.path(data_dir, panel, "source_metrics.tsv"))
  source <- source[source$target %in% d$target, ]
  original <- d
  original$method <- original$source_method
  check_rows(original, source, c("method", "target"), c("Pearson", "status"))
  # Sector area encodes |r|; its color retains the sign. Missing values stay missing.
  sector <- abs(d$Pearson)
  symbol <- ifelse(d$status == "valid", "sector", ifelse(d$status == "unavailable", "/", "x"))
  same_numbers(sector, d$sector_fraction)
  stopifnot(identical(symbol, d$symbol))
  d$sector_fraction <- sector
  d$symbol <- symbol
  d
}
verify_analysis <- function(analysis, data_dir) {
  d <- do.call(rbind, analysis[LETTERS[1:4]])
  stopifnot(nrow(d) == 469L, sum(d$status == "valid") == 332L)
  reference <- read_table(file.path(data_dir, "Fig2_all_panels_correlations.tsv"))
  check_rows(d, reference, c("dataset", "method", "target"), names(d))
  points <- read_table(file.path(data_dir, "E", "mean_SD_points.tsv"))
  check_rows(d[d$status == "valid", ], points, c("dataset", "method", "target"), "Pearson")
  cat("Fig2: 469 correlation records and 28 recalculated mean/SD summaries verified.\n")
}
