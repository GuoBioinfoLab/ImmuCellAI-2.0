args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/age_associated_remodeling/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2"))

fraction <- read_fraction_matrix(file.path(output_dir, "deconvolution", "state_fraction.txt"))
source(file.path(repo_root, "analysis", "common", "age_groups.R"))
missing <- setdiff(unlist(figure5_major_cells), colnames(fraction))
if (length(missing)) stop("Missing Figure 5 states: ", paste(missing, collapse = ", "))
major <- vapply(figure5_major_cells, function(cells) {
  rowSums(fraction[, cells, drop = FALSE])
}, numeric(nrow(fraction)))
rownames(major) <- rownames(fraction)
meta <- read_table_auto(file.path(output_dir, "healthy_age_aligned_metadata.tsv"))
age_group <- figure5_age_group
levels_age <- figure5_age_levels
dat <- data.frame(sample = rownames(major), major, check.names = FALSE) |>
  dplyr::left_join(meta, by = "sample") |>
  dplyr::mutate(AgeGroup = factor(age_group(age), levels = levels_age)) |>
  dplyr::filter(!is.na(AgeGroup)) |>
  tidyr::pivot_longer(dplyr::all_of(colnames(major)), names_to = "Lineage", values_to = "Fraction") |>
  dplyr::group_by(AgeGroup, Lineage) |>
  dplyr::summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop") |>
  dplyr::group_by(AgeGroup) |>
  dplyr::mutate(Fraction = Fraction / sum(Fraction)) |>
  dplyr::ungroup() |>
  dplyr::mutate(Lineage = factor(Lineage, levels = names(figure5_major_cells)))
data.table::fwrite(dat, file.path(output_dir, "age_major_lineage_composition.tsv"), sep = "\t")
p <- ggplot2::ggplot(dat, ggplot2::aes(AgeGroup, Fraction, fill = Lineage)) +
  ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.2) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_fill_manual(
    breaks = names(figure5_major_cells),
    values = c(Tcell = "#E15759", Bcell = "#4E79A7", ILC = "#F28E2B", NK = "#76B7B2",
               Granulocyte = "#59A14F", Monocyte = "#EDC948", Macrophage = "#B07AA1", DC = "#FF9DA7")
  ) +
  ggplot2::theme_bw() + ggplot2::labs(x = "Age group", y = "Mean immune composition", fill = NULL) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_plot_pair(p, file.path(output_dir, "Age_major_lineage_stacked_barplot"), 6.5, 5)
