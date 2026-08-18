args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/age_associated_remodeling/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2"))

major <- read_fraction_matrix(file.path(output_dir, "deconvolution", "major_lineage_fraction.txt"))
meta <- read_table_auto(file.path(output_dir, "healthy_age_aligned_metadata.tsv"))
age_group <- function(age) dplyr::case_when(
  age >= 0 & age <= 1 ~ "0-1", age >= 10 & age <= 20 ~ "10-20",
  age > 20 & age <= 30 ~ "20-30", age > 30 & age <= 50 ~ "30-50",
  age > 50 & age <= 70 ~ "50-70", age > 70 & age <= 90 ~ "70-90",
  age > 90 ~ "90+", TRUE ~ NA_character_
)
levels_age <- c("0-1", "10-20", "20-30", "30-50", "50-70", "70-90", "90+")
dat <- data.frame(sample = rownames(major), major, check.names = FALSE) |>
  dplyr::left_join(meta, by = "sample") |>
  dplyr::mutate(AgeGroup = factor(age_group(age), levels = levels_age)) |>
  dplyr::filter(!is.na(AgeGroup)) |>
  tidyr::pivot_longer(dplyr::all_of(colnames(major)), names_to = "Lineage", values_to = "Fraction") |>
  dplyr::group_by(AgeGroup, Lineage) |>
  dplyr::summarise(Fraction = mean(Fraction, na.rm = TRUE), .groups = "drop") |>
  dplyr::group_by(AgeGroup) |>
  dplyr::mutate(Fraction = Fraction / sum(Fraction)) |>
  dplyr::ungroup()
data.table::fwrite(dat, file.path(output_dir, "age_major_lineage_composition.tsv"), sep = "\t")
p <- ggplot2::ggplot(dat, ggplot2::aes(AgeGroup, Fraction, fill = Lineage)) +
  ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.2) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::theme_bw() + ggplot2::labs(x = "Age group", y = "Mean immune composition", fill = NULL) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_plot_pair(p, file.path(output_dir, "Age_major_lineage_stacked_barplot"), 6.5, 5)
