# Aggregate 53 cell types per sample, then take means within each age group.
mapping <- read_tsv(file.path(data_dir, "D/major_cell_mapping.tsv"))
stopifnot(nrow(mapping) == 53L, !anyDuplicated(mapping$source_cell),
          setequal(mapping$source_cell, colnames(fraction)))
major_names <- unique(mapping$major_cell)
major <- sapply(major_names, function(cell) rowSums(fraction[, mapping$source_cell[mapping$major_cell == cell], drop = FALSE]))
stored_major <- read_tsv(file.path(data_dir, "D/sample_major_fractions.tsv"))
stopifnot(identical(stored_major$Sample.ID, meta$Sample.ID), identical(stored_major$group, meta$group))
compare("D_sample_major_fractions", unname(as.matrix(stored_major[, major_names])), unname(major))
means <- do.call(rbind, lapply(age_levels, function(g) colMeans(major[meta$group == g, , drop = FALSE])))
rownames(means) <- age_levels
stored_means <- read_tsv(file.path(data_dir, "D/age_group_major_means.tsv"))
stopifnot(identical(stored_means$group, age_levels))
compare("D_group_means", unname(as.matrix(stored_means[, major_names])), unname(means))
bars <- read_tsv(file.path(data_dir, "D/plot_data.tsv"))
stopifnot(nrow(bars) == 48L, !anyDuplicated(paste(bars$group, bars$cell_type)))
compare("D_barplot_values", bars$proportion, unname(means[cbind(match(bars$group, age_levels), match(bars$cell_type, major_names))]))
compare("D_barplot_totals", unname(rowSums(means)), rep(1, 6))
analysis$D <- bars
