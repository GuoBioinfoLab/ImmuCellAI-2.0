# Shared input checks and the original age-group/IQR display rules.
age_levels <- c("age0-1", "age10-20", "age20-30", "age30-50", "age50-70", "age70+")
panel_cells <- list(
  A = c("CD4Tn", "CD8Tn", "M0", "gdT"),
  B = c("cNK", "CD8Temra", "ncMo", "Eosinophil"),
  C = c("CD4Tem", "CD8Tem", "CD4Tcm", "Basophil", "cMo", "Breg", "Bn", "MZB", "ILC1")
)
read_tsv <- function(path, matrix = FALSE) {
  if (matrix) return(as.matrix(read.delim(path, row.names = 1, check.names = FALSE,
                                        stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")))
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}
numeric_errors <- list()
compare <- function(label, observed, expected, tolerance = 1e-12) {
  stopifnot(identical(dim(observed), dim(expected)), length(observed) == length(expected),
            all(is.finite(observed)), all(is.finite(expected)))
  difference <- max(abs(observed - expected), 0)
  if (difference > tolerance) stop(label, " differs: ", difference)
  numeric_errors[[label]] <<- difference
}
key <- function(x) paste(x$Sample.ID, x$cell_type, sep = "|")
verify_inputs <- function(data_dir) {
  fraction <- read_tsv(file.path(data_dir, "_shared/state_fraction_sample_by_celltype.tsv"), TRUE)
  meta <- read_tsv(file.path(data_dir, "_shared/sample_metadata.tsv"))
  stopifnot(identical(dim(fraction), c(327L, 53L)), nrow(meta) == 327L,
            !anyDuplicated(rownames(fraction)), !anyDuplicated(meta$Sample.ID),
            identical(rownames(fraction), meta$Sample.ID), all(fraction >= 0))
  compare("sample_fraction_totals", unname(rowSums(fraction)), rep(1, 327))
  expected_group <- rep(NA_character_, nrow(meta))
  age <- meta$Age
  expected_group[age >= 0 & age < 2] <- "age0-1"
  expected_group[age > 10 & age <= 20] <- "age10-20"
  expected_group[age > 20 & age <= 30] <- "age20-30"
  expected_group[age > 30 & age < 50] <- "age30-50"
  expected_group[age >= 50 & age < 70] <- "age50-70"
  expected_group[age >= 70] <- "age70+"
  stopifnot(identical(meta$group, expected_group))
  counts <- as.integer(table(factor(meta$group, levels = age_levels)))
  stopifnot(identical(counts, c(7L, 43L, 124L, 85L, 60L, 8L)))
  long <- read_tsv(file.path(data_dir, "_shared/all_53_cells_long.tsv"))
  flags <- read_tsv(file.path(data_dir, "_shared/all_53_cells_with_outlier_flags.tsv"))
  key <- function(x) paste(x$Sample.ID, x$cell_type, sep = "|")
  stopifnot(nrow(long) == 327L * 53L, !anyDuplicated(key(long)), identical(key(long), key(flags)))
  source_cells <- long$cell_type
  source_cells[source_cells == "Bn"] <- "Bnaive"
  source_cells[source_cells == "moDC"] <- "monoDC"
  original_values <- fraction[cbind(match(long$Sample.ID, rownames(fraction)), match(source_cells, colnames(fraction)))]
  compare("all_long_values", long$value, unname(original_values))
  compare("flagged_values", flags$value, long$value)
  stopifnot(identical(long$group, meta$group[match(long$Sample.ID, meta$Sample.ID)]))
  recomputed <- logical(nrow(long))
  groups <- split(seq_len(nrow(long)), interaction(long$cell_type, long$group, drop = TRUE))
  for (idx in groups) {
    value <- long$value[idx]
    q <- quantile(value, c(.25, .75), type = 7)
    iqr <- q[2] - q[1]
    recomputed[idx] <- if (iqr <= 0) FALSE else value < q[1] - 1.5 * iqr | value > q[2] + 1.5 * iqr
    compare(paste0("quartiles_", long$cell_type[idx[1]], "_", long$group[idx[1]]),
            c(flags$plot_q1[idx[1]], flags$plot_q3[idx[1]]), unname(q))
  }
  stopifnot(identical(recomputed, flags$plot_outlier))
  old_filtered <- read_tsv(file.path(data_dir, "_shared/original_filtered_long.tsv"))
  stopifnot(identical(key(old_filtered), key(flags[!recomputed | flags$cell_type == "CD8Tem", ])))
  list(fraction = fraction, meta = meta, flags = flags)
}
analyze_box_panel <- function(panel, flags, data_dir) {
  raw <- read_tsv(file.path(data_dir, panel, "all_values_with_display_flags.tsv"))
  plotted <- read_tsv(file.path(data_dir, panel, "plot_data.tsv"))
  selected <- flags[flags$cell_type %in% panel_cells[[panel]], ]
  filter_cell <- selected$cell_type == "CD8Temra" |
    (panel == "C" & !selected$cell_type %in% c("CD4Tem", "CD8Tem"))
  keep <- !(filter_cell & selected$plot_outlier)
  stopifnot(identical(key(raw), key(selected)), identical(raw$include_in_plot, keep),
            identical(key(plotted), key(selected[keep, ])))
  compare(paste0(panel, "_plotted_values"), plotted$value, selected$value[keep])
  order <- read_tsv(file.path(data_dir, panel, "cell_order.tsv"))
  stopifnot(identical(order$cell_type, panel_cells[[panel]]))
  panel_counts[[panel]] <<- list(n_cells = length(panel_cells[[panel]]), raw = nrow(raw),
                                plotted = nrow(plotted), display_only_excluded = sum(!keep))
  plotted
}
