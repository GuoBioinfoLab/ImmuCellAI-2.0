# Mean and sample SD across each tool's evaluable displayed targets.
correlations <- do.call(rbind, analysis[LETTERS[1:4]])
rows <- list()
for (dataset in datasets) for (method in methods) {
  d <- correlations[correlations$dataset == dataset & correlations$method == method, ]
  r <- d$Pearson[d$status == "valid"]
  stopifnot(length(r) >= 2L)
  avg <- mean(r)
  std <- sd(r)
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = dataset, scope = "all_valid", method = method,
    n = length(r), n_candidate = nrow(d), mean_Pearson = avg, SD_Pearson = std,
    lower = avg - std, upper = avg + std)
}
analysis$E <- do.call(rbind, rows)
expected <- read_table(file.path(data_dir, "E", "mean_SD_summary.tsv"))
check_rows(analysis$E, expected, c("dataset", "method"), names(analysis$E))
