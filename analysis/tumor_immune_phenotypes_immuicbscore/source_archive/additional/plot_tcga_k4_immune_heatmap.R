library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- file.path(fig4_dir, "cluster_TCGA")
state_file <- file.path(fig4_dir, "TCGA_ImmuCellAI2_state_fraction.txt")
class_file <- file.path(cluster_dir, "cluster_TCGA.k=4.consensusClass.csv")

out_sample_pdf <- file.path(cluster_dir, "TCGA_ImmuCellAI2_K4_cell_fraction_sample_heatmap.pdf")
out_sample_png <- file.path(cluster_dir, "TCGA_ImmuCellAI2_K4_cell_fraction_sample_heatmap.png")
out_mean_pdf <- file.path(cluster_dir, "TCGA_ImmuCellAI2_K4_cluster_mean_cell_fraction_heatmap.pdf")
out_mean_png <- file.path(cluster_dir, "TCGA_ImmuCellAI2_K4_cluster_mean_cell_fraction_heatmap.png")
out_feature_table <- file.path(cluster_dir, "TCGA_ImmuCellAI2_K4_cluster_celltype_features.txt")
out_aligned_matrix <- file.path(cluster_dir, "TCGA_ImmuCellAI2_K4_aligned_state_fraction.txt")

dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading state fraction matrix...")
state_dt <- fread(state_file, check.names = FALSE)
sample_ids <- state_dt[[1L]]
state_dt[[1L]] <- NULL
state_mat <- as.matrix(state_dt)
storage.mode(state_mat) <- "double"
rownames(state_mat) <- sample_ids
rm(state_dt)
gc()

message("Reading k=4 consensus class...")
class_dt <- fread(class_file, header = FALSE, col.names = c("Sample", "Cluster"))
class_dt[, Cluster := factor(paste0("C", Cluster), levels = paste0("C", 1:4))]

common_samples <- intersect(class_dt$Sample, rownames(state_mat))
if (!length(common_samples)) stop("No overlapping samples between state matrix and k=4 class file.")

state_mat <- state_mat[common_samples, , drop = FALSE]
class_dt <- class_dt[match(common_samples, Sample)]

complete_rows <- rowSums(is.na(state_mat)) == 0
if (any(!complete_rows)) {
  message("Removing samples with NA predictions: ", sum(!complete_rows))
}
state_mat <- state_mat[complete_rows, , drop = FALSE]
class_dt <- class_dt[complete_rows]

fwrite(
  data.table(Sample = rownames(state_mat), state_mat, check.names = FALSE),
  out_aligned_matrix,
  sep = "\t",
  quote = FALSE
)

message("Preparing z-score matrix...")
# Heatmap rows are cell types, columns are samples.
heat_mat_raw <- t(state_mat)
heat_mat_z <- t(scale(t(heat_mat_raw)))
heat_mat_z[!is.finite(heat_mat_z)] <- 0
heat_mat_z <- pmax(pmin(heat_mat_z, 3), -3)

cluster_colors <- c(
  C1 = "#E64B35",
  C2 = "#4DBBD5",
  C3 = "#00A087",
  C4 = "#3C5488"
)
heat_col <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))

cluster_order <- levels(class_dt$Cluster)
column_order <- unlist(lapply(cluster_order, function(cl) {
  idx <- which(class_dt$Cluster == cl)
  if (length(idx) <= 2L) return(idx)
  # Order samples within each consensus cluster by their cell-fraction pattern.
  idx[hclust(dist(t(heat_mat_z[, idx, drop = FALSE])), method = "average")$order]
}), use.names = FALSE)

ha <- HeatmapAnnotation(
  Cluster = class_dt$Cluster[column_order],
  col = list(Cluster = cluster_colors),
  annotation_name_gp = gpar(fontsize = 9),
  show_annotation_name = TRUE
)

sample_heatmap <- Heatmap(
  heat_mat_z[, column_order, drop = FALSE],
  name = "Cell fraction\nrow z-score",
  col = heat_col,
  top_annotation = ha,
  column_split = class_dt$Cluster[column_order],
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title = "Immune cell types",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 2,
  border = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

message("Writing sample-level heatmap...")
pdf(out_sample_pdf, width = 13, height = 9)
draw(sample_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
png(out_sample_png, width = 3900, height = 2700, res = 300)
draw(sample_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Computing cluster means and features...")
cluster_mean <- sapply(cluster_order, function(cl) {
  colMeans(state_mat[class_dt$Cluster == cl, , drop = FALSE], na.rm = TRUE)
})
cluster_mean_z <- t(scale(t(cluster_mean)))
cluster_mean_z[!is.finite(cluster_mean_z)] <- 0
cluster_mean_z <- pmax(pmin(cluster_mean_z, 3), -3)

feature_dt <- rbindlist(lapply(cluster_order, function(cl) {
  target <- cluster_mean[, cl]
  others <- rowMeans(cluster_mean[, setdiff(cluster_order, cl), drop = FALSE])
  z_target <- cluster_mean_z[, cl]
  data.table(
    Cluster = cl,
    CellType = rownames(cluster_mean),
    MeanFraction = target,
    MeanOtherClusters = others,
    Difference = target - others,
    ZScoreAcrossClusters = z_target
  )[order(-ZScoreAcrossClusters, -Difference)]
}), use.names = TRUE)
fwrite(feature_dt, out_feature_table, sep = "\t", quote = FALSE)

mean_top <- HeatmapAnnotation(
  Cluster = factor(colnames(cluster_mean_z), levels = cluster_order),
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE
)

mean_heatmap <- Heatmap(
  cluster_mean_z,
  name = "Cluster mean\nrow z-score",
  col = heat_col,
  top_annotation = mean_top,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title = "Immune cell types",
  column_title = "K = 4 consensus clusters",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", cluster_mean[i, j]), x, y, gp = gpar(fontsize = 5.5))
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

message("Writing cluster-mean heatmap...")
pdf(out_mean_pdf, width = 6.5, height = 9)
draw(mean_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
png(out_mean_png, width = 1950, height = 2700, res = 300)
draw(mean_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Top enriched cell types by cluster:")
for (cl in cluster_order) {
  top <- feature_dt[Cluster == cl][1:10, .(CellType, MeanFraction, Difference, ZScoreAcrossClusters)]
  message("\n", cl)
  print(top)
}

message("Saved files:")
message(out_sample_pdf)
message(out_sample_png)
message(out_mean_pdf)
message(out_mean_png)
message(out_feature_table)
