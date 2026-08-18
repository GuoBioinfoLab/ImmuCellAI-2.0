args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/tumor_immune_phenotypes_immuicbscore/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "ComplexHeatmap", "circlize", "survival", "survminer", "ggplot2"))

fraction <- read_fraction_matrix(tcga_fraction_file)
cluster <- read_table_auto(tcga_cluster_file)
if (!all(c("sample", "cluster") %in% names(cluster))) stop("Cluster file requires sample and cluster columns.")
cluster$cluster <- factor(cluster$cluster, levels = paste0("IP", 1:4))
cluster <- cluster[match(rownames(fraction), cluster$sample), , drop = FALSE]
if (anyNA(cluster$cluster)) stop("Some fraction samples are absent from the K=4 class file.")

set.seed(123)
order_samples <- unlist(lapply(levels(cluster$cluster), function(z) {
  sample(cluster$sample[cluster$cluster == z])
}))
heat <- zscore_rows(fraction[order_samples, , drop = FALSE], clip = 1)

cluster_colors <- c(IP1 = "#E15759", IP2 = "#4E79A7", IP3 = "#F28E2B", IP4 = "#76B7B2")
column_annotation <- ComplexHeatmap::HeatmapAnnotation(
  ImmunePhenotype = factor(cluster$cluster[match(order_samples, cluster$sample)], levels = names(cluster_colors)),
  col = list(ImmunePhenotype = cluster_colors),
  show_annotation_name = FALSE
)
ht <- ComplexHeatmap::Heatmap(
  heat,
  name = "Row z-score",
  col = circlize::colorRamp2(c(-1, 0, 1), c("#3341A3", "black", "#E0DA54")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = factor(cluster$cluster[match(order_samples, cluster$sample)], levels = names(cluster_colors)),
  top_annotation = column_annotation,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  border = TRUE
)
pdf(file.path(output_dir, "TCGA_four_immune_phenotypes_heatmap.pdf"), width = 8, height = 10)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

surv <- read_table_auto(tcga_survival_file)
required <- c("sample", "time_months", "status")
if (!all(required %in% names(surv))) {
  stop("Survival table must contain sample, time_months, and status (0=censored, 1=event).")
}
surv$case_id <- substr(as.character(surv$sample), 1L, 12L)
cluster$case_id <- substr(as.character(cluster$sample), 1L, 12L)
dat <- merge(cluster[, c("case_id", "cluster")], surv, by = "case_id")
dat <- dat[is.finite(dat$time_months) & dat$status %in% c(0, 1), , drop = FALSE]
fit <- survival::survfit(survival::Surv(time_months, status) ~ cluster, data = dat)
p <- survminer::ggsurvplot(
  fit, data = dat, pval = TRUE, risk.table = TRUE,
  palette = cluster_colors, xlab = "Months", ylab = "Overall survival",
  legend.title = "Immune phenotype"
)
pdf(file.path(output_dir, "TCGA_four_immune_phenotypes_overall_survival.pdf"), width = 6.5, height = 6.2)
print(p)
dev.off()
data.table::fwrite(dat, file.path(output_dir, "TCGA_survival_analysis_table.tsv"), sep = "\t")
