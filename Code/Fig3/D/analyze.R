# Sample SD across targets, not standard error or a confidence interval.
records <- do.call(rbind, analysis[LETTERS[1:3]])
rows <- list()
for (dataset in groups) for (method in methods) {
  d <- records[records$Dataset == dataset & records$method == method, ]
  valid <- d[is.finite(d$r), ]
  stopifnot(nrow(valid) >= 2L)
  avg <- mean(valid$r)
  std <- sd(valid$r)
  rows[[length(rows) + 1L]] <- data.frame(
    Dataset = dataset, source_panel = LETTERS[match(dataset, groups)], method = method,
    n_valid = nrow(valid), n_candidates = nrow(d), mean = avg, SD = std,
    lower = avg - std, upper = avg + std,
    included_targets = paste(valid$cellType, collapse = " | "))
}
analysis$D <- do.call(rbind, rows)
expected <- read_csv(file.path(data_dir, "D", "mean_SD_summary.csv"))
key <- function(d) paste(d$Dataset, d$method, sep = "\r")
stopifnot(nrow(expected) == 21L, !anyDuplicated(key(expected)),
          setequal(key(analysis$D), key(expected)))
expected <- expected[match(key(analysis$D), key(expected)), ]
for (field in c("n_valid", "n_candidates", "mean", "SD", "lower", "upper")) {
  same_numbers(analysis$D[[field]], as.numeric(expected[[field]]))
}
stopifnot(identical(analysis$D$included_targets, expected$included_targets))
