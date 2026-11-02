# Recompute group means, sample-SD z-scores, clipping and the archived column order.
all_means <- do.call(rbind, lapply(age_levels, function(g) colMeans(fraction[meta$group == g, , drop = FALSE])))
rownames(all_means) <- age_levels
hm_means <- read_tsv(file.path(data_dir, "E/age_group_mean_fraction.tsv"), TRUE)
compare("E_group_means", unname(hm_means), unname(all_means[, colnames(hm_means)]))
scaled <- scale(hm_means)
scaled[is.na(scaled)] <- 0
scaled <- pmax(-1, pmin(1, scaled))
# pmin/pmax keep dimensions only when the matrix is the first argument.
dim(scaled) <- dim(hm_means)
dimnames(scaled) <- dimnames(hm_means)
peak <- apply(scaled, 2, which.max)
category <- ifelse(peak == 1, 1, ifelse(peak %in% 2:4, 2, 3))
peak_value <- scaled[cbind(peak, seq_along(peak))]
ordered <- colnames(scaled)[order(category, peak, -peak_value, seq_along(peak))]
first <- intersect(c("CD4Tn", "CD8Tn", "M0", "gdT"), ordered)
ordered <- c(first, setdiff(ordered, first))
swap <- match(c("FOB", "CD8Temra"), ordered)
ordered[swap] <- rev(ordered[swap])
stored_order <- read_tsv(file.path(data_dir, "E/cell_order.tsv"))
stopifnot(identical(ordered, stored_order$cell_type), identical(stored_order$final_rank, seq_len(48)))
hm <- read_tsv(file.path(data_dir, "E/plot_matrix_zscore_clipped.tsv"), TRUE)
labels <- ordered
labels[labels == "Bnaive"] <- "Bn"
labels[labels == "monoDC"] <- "moDC"
labels[labels == "gdT"] <- paste0(intToUtf8(c(947, 948)), "T")
stopifnot(identical(dim(hm), c(6L, 48L)), identical(rownames(hm), age_levels), identical(colnames(hm), labels))
compare("E_clipped_zscores", unname(hm), unname(scaled[, ordered]))
analysis$E <- hm
