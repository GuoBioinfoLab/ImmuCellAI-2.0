args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/tumor_immune_phenotypes_immuicbscore/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "ranger", "ggplot2", "dplyr", "tidyr"))

icb_fraction <- read_fraction_matrix(icb_fraction_file)
tcga_fraction <- read_fraction_matrix(tcga_fraction_file)
common_cells <- intersect(colnames(icb_fraction), colnames(tcga_fraction))
if (length(common_cells) != 53L) warning("Using ", length(common_cells), " shared immune states.")
clinical <- read_table_auto(icb_clinical_file)
if (!all(c("Run", "SRA_study", "ResponseBinary") %in% names(clinical))) {
  stop("ICB clinical table requires Run, SRA_study, and ResponseBinary.")
}
clinical <- clinical[clinical$SRA_study %in% selected_icb_studies & clinical$ResponseBinary %in% c(0, 1), ]
idx <- match(clinical$Run, rownames(icb_fraction))
keep <- !is.na(idx)
clinical <- clinical[keep, ]
x <- clr_transform(icb_fraction[idx[keep], common_cells, drop = FALSE])
center <- colMeans(x)
scale_value <- apply(x, 2, sd)
scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
x <- sweep(sweep(x, 2, center, "-"), 2, scale_value, "/")
y <- as.integer(clinical$ResponseBinary)
weights <- c("0" = 0.5 / mean(y == 0), "1" = 0.5 / mean(y == 1))
set.seed(9001)
model <- ranger::ranger(
  x = as.data.frame(x), y = factor(y, levels = c(0, 1)),
  probability = TRUE, num.trees = 1200,
  mtry = max(1L, floor(sqrt(ncol(x)))), min.node.size = 8,
  class.weights = weights, num.threads = n_cores, seed = 9001
)

tcga_x <- clr_transform(tcga_fraction[, common_cells, drop = FALSE])
tcga_x <- sweep(sweep(tcga_x, 2, center, "-"), 2, scale_value, "/")
score <- predict(model, data = as.data.frame(tcga_x), num.threads = n_cores)$predictions[, "1"]
cluster <- read_table_auto(tcga_cluster_file)
score_table <- data.frame(sample = rownames(tcga_fraction), ImmuICBscore = score)
score_table <- merge(score_table, cluster, by = "sample")
score_table$cluster <- factor(score_table$cluster, levels = paste0("IP", 1:4))
score_table$quartile <- cut(
  score_table$ImmuICBscore,
  quantile(score_table$ImmuICBscore, seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE, labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")
)
data.table::fwrite(score_table, file.path(output_dir, "TCGA_ImmuICBscore_and_phenotype.tsv"), sep = "\t")

cols <- c(IP1 = "#E15759", IP2 = "#4E79A7", IP3 = "#F28E2B", IP4 = "#76B7B2")
p1 <- ggplot2::ggplot(score_table, ggplot2::aes(cluster, ImmuICBscore, fill = cluster)) +
  ggplot2::geom_violin(trim = FALSE) + ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA) +
  ggplot2::scale_fill_manual(values = cols) + ggplot2::theme_bw() +
  ggplot2::labs(x = "TCGA immune phenotype", y = "ImmuICBscore") +
  ggplot2::theme(legend.position = "none")

composition <- score_table |>
  dplyr::count(quartile, cluster) |>
  dplyr::group_by(quartile) |>
  dplyr::mutate(proportion = n / sum(n)) |>
  dplyr::ungroup()
p2 <- ggplot2::ggplot(composition, ggplot2::aes(quartile, proportion, fill = cluster)) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::theme_bw() + ggplot2::labs(x = "ImmuICBscore quartile", y = "Cluster proportion", fill = "Phenotype")
save_plot_pair(p1, file.path(output_dir, "TCGA_ImmuICBscore_by_phenotype"), 4.6, 4.2)
save_plot_pair(p2, file.path(output_dir, "TCGA_phenotype_composition_by_score_quartile"), 5.5, 4.2)
