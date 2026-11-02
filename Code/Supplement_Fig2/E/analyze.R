
d <- read_tsv(file.path(data_dir, "E/mean_SD_plot_data.tsv"))
stopifnot(nrow(d) == 28L, !anyDuplicated(table_key(d, c("dataset", "method"))))
check_table(d, read_tsv(file.path(data_dir, "_shared/mean_SD_plot_data.tsv")), c("dataset", "method"))
for (i in seq_len(nrow(d))) {
 row <- d[i, ]; z <- context$shown[context$shown$dataset == row$dataset & context$shown$method == row$method, ]
 r <- z$Pearson[z$status == "valid"]
 stopifnot(all(is.finite(r)), all(is.na(z$Pearson[z$status != "valid"])), length(r) == row$n_valid, nrow(z) == row$n_targets)
 check_num(c(row$mean_Pearson, row$SD_Pearson), c(mean(r), sd(r)), "Mean and sample SD")
}
analysis$E <- d
