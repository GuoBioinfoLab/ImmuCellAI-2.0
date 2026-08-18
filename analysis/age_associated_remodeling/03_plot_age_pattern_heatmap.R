args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/age_associated_remodeling/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "dplyr", "ComplexHeatmap", "circlize"))

fraction <- read_fraction_matrix(file.path(output_dir, "deconvolution", "state_fraction.txt"))
meta <- read_table_auto(file.path(output_dir, "healthy_age_aligned_metadata.tsv"))
age_group <- function(age) dplyr::case_when(
  age >= 0 & age <= 1 ~ "0-1", age >= 10 & age <= 20 ~ "10-20",
  age > 20 & age <= 30 ~ "20-30", age > 30 & age <= 50 ~ "30-50",
  age > 50 & age <= 70 ~ "50-70", age > 70 & age <= 90 ~ "70-90",
  age > 90 ~ "90+", TRUE ~ NA_character_
)
levels_age <- c("0-1", "10-20", "20-30", "30-50", "50-70", "70-90", "90+")
dat <- data.frame(sample = rownames(fraction), fraction, check.names = FALSE) |>
  dplyr::left_join(meta, by = "sample") |>
  dplyr::mutate(AgeGroup = factor(age_group(age), levels = levels_age)) |>
  dplyr::filter(!is.na(AgeGroup))
mean_by_group <- dat |>
  dplyr::group_by(AgeGroup) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(colnames(fraction)), mean, na.rm = TRUE), .groups = "drop")
mat <- as.matrix(mean_by_group[, colnames(fraction), drop = FALSE])
rownames(mat) <- as.character(mean_by_group$AgeGroup)
z <- scale(mat)
z[!is.finite(z)] <- 0
z[z > 1] <- 1
z[z < -1] <- -1

fixed_first <- intersect(c("CD4Tn", "CD8Tn", "M0", "gdT"), colnames(z))
remaining <- setdiff(colnames(z), fixed_first)
peak_group <- apply(z[, remaining, drop = FALSE], 2, which.max)
remaining <- remaining[order(peak_group, colnames(z)[match(remaining, colnames(z))])]
cell_order <- c(fixed_first, remaining)
z <- z[, cell_order, drop = FALSE]
data.table::fwrite(data.frame(AgeGroup = rownames(z), z, check.names = FALSE), file.path(output_dir, "age_group_mean_scaled.tsv"), sep = "\t")

ht <- ComplexHeatmap::Heatmap(
  z, name = "Scaled mean",
  col = circlize::colorRamp2(c(-1, 0, 1), c("#3341A3", "black", "#E0DA54")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  column_names_rot = 90,
  heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1))
)
pdf(file.path(output_dir, "Age_associated_immune_pattern_heatmap.pdf"), width = 12, height = 4)
ComplexHeatmap::draw(ht)
dev.off()
