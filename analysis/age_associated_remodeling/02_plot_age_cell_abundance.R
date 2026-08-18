args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/age_associated_remodeling/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2", "patchwork"))

fraction <- read_fraction_matrix(file.path(output_dir, "deconvolution", "state_fraction.txt"))
meta <- read_table_auto(file.path(output_dir, "healthy_age_aligned_metadata.tsv"))

age_group <- function(age) {
  dplyr::case_when(
    age >= 0 & age <= 1 ~ "0-1",
    age >= 10 & age <= 20 ~ "10-20",
    age > 20 & age <= 30 ~ "20-30",
    age > 30 & age <= 50 ~ "30-50",
    age > 50 & age <= 70 ~ "50-70",
    age > 70 & age <= 90 ~ "70-90",
    age > 90 ~ "90+",
    TRUE ~ NA_character_
  )
}
group_levels <- c("0-1", "10-20", "20-30", "30-50", "50-70", "70-90", "90+")
dat <- data.frame(sample = rownames(fraction), fraction, check.names = FALSE) |>
  dplyr::left_join(meta, by = "sample") |>
  dplyr::mutate(AgeGroup = factor(age_group(age), levels = group_levels)) |>
  dplyr::filter(!is.na(AgeGroup))

panel_a <- c("CD4Tn", "CD8Tn", "M0", "gdT")
panel_b <- c("cNK", "CD8Temra", "ncMo", "Eosinophil")
panel_c <- c("CD4Tem", "CD8Tem", "CD4Tcm", "Basophil", "cMo", "Breg", "Bnaive", "MZB", "ILC1")
all_cells <- unique(c(panel_a, panel_b, panel_c))
missing <- setdiff(all_cells, names(dat))
if (length(missing)) stop("Missing selected cells: ", paste(missing, collapse = ", "))

long <- tidyr::pivot_longer(dat, dplyr::all_of(all_cells), names_to = "CellType", values_to = "Fraction")
data.table::fwrite(long, file.path(output_dir, "age_selected_cells_long.tsv"), sep = "\t")
palette <- c("#C95D5D", "#5B7FA8", "#E69F3A", "#83A66A", "#62A6A0", "#7B6AA6", "#E0C45C")
names(palette) <- group_levels

trim_for_display <- function(x) {
  q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
  upper <- q[2] + 1.5 * diff(q)
  pmin(x, upper)
}

make_one <- function(cell, trim = FALSE, upper = NULL) {
  x <- long[long$CellType == cell, , drop = FALSE]
  x$DisplayFraction <- if (trim) trim_for_display(x$Fraction) else x$Fraction
  if (!is.null(upper)) x$DisplayFraction <- pmin(x$DisplayFraction, upper)
  means <- x |>
    dplyr::group_by(AgeGroup) |>
    dplyr::summarise(DisplayFraction = mean(DisplayFraction, na.rm = TRUE), .groups = "drop")
  ggplot2::ggplot(x, ggplot2::aes(AgeGroup, DisplayFraction, fill = AgeGroup)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.72, linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.18, size = 0.7, shape = 21, fill = "white", alpha = 0.7) +
    ggplot2::geom_line(data = means, ggplot2::aes(group = 1), color = "#4E79A7", linewidth = 0.7) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(title = cell, x = "Age group", y = "Cell fraction") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), aspect.ratio = 1)
}

plots_a <- lapply(panel_a, function(z) make_one(z, upper = if (z == "M0") 0.018 else NULL))
plots_b <- lapply(panel_b, function(z) make_one(z, trim = z == "CD8Temra"))
plots_c <- lapply(panel_c, function(z) make_one(
  z, trim = !z %in% c("CD4Tem", "CD8Tem"), upper = if (z == "CD8Tem") 0.0038 else NULL
))
patchwork::wrap_plots(plots_a, nrow = 1) |> save_plot_pair(file.path(output_dir, "Age_abundance_panel_A"), 12, 3.2)
patchwork::wrap_plots(plots_b, nrow = 1) |> save_plot_pair(file.path(output_dir, "Age_abundance_panel_B"), 12, 3.2)
patchwork::wrap_plots(plots_c, ncol = 3) |> save_plot_pair(file.path(output_dir, "Age_abundance_panel_C"), 9, 9)
