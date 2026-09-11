options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
  library(grid)
})

fig5_dir <- "<LOCAL_R_ROOT>/Fig5"
out_dir <- file.path(fig5_dir, "ImmuCellAI2_age_group_abundance")

fraction_file <- file.path(out_dir, "ImmuCellAI2_age_state_fraction_sample_by_celltype.txt")
sample_info_file <- file.path(out_dir, "ImmuCellAI2_age_combined_sample_info_with_groups.txt")

age_levels <- c("age0-1", "age10-20", "age20-30", "age30-50", "age50-70", "age70+")

message("Reading ImmuCellAI2 cell fraction matrix...")
fraction <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
sample_col <- names(fraction)[1]
rownames(fraction) <- fraction[[sample_col]]
fraction[[sample_col]] <- NULL
fraction[] <- lapply(fraction, function(x) as.numeric(as.character(x)))
fraction <- as.data.frame(fraction, check.names = FALSE)

message("Reading age group metadata...")
sample_info <- fread(sample_info_file, data.table = FALSE, check.names = FALSE)
sample_info <- sample_info %>%
  mutate(group = factor(group, levels = age_levels)) %>%
  filter(Sample.ID %in% rownames(fraction), !is.na(group))

fraction <- fraction[sample_info$Sample.ID, , drop = FALSE]
fraction$group <- sample_info$group

message("Calculating group means...")
grouped_mean <- fraction %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), function(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

rownames(grouped_mean) <- as.character(grouped_mean$group)
grouped_mean$group <- NULL
grouped_mean <- grouped_mean[age_levels[age_levels %in% rownames(grouped_mean)], , drop = FALSE]

heatmap_cells <- c(
  "cMo", "M0", "monoDC", "CD4Tn", "CD8Tn", "Mast cell", "ILC2", "Basophil",
  "Breg", "Neutrophil", "CD4Tem", "CD8Tem", "FOB", "BGC", "gdT", "ILC1",
  "MAIT", "MBC", "MZB", "NKT", "Th1", "Th2", "Tfh", "Bnaive", "Treg", "ncMo",
  "pDC", "Tr1", "CD4Temra", "CD8Temra", "cDC1", "cDC2", "Bex", "Th1/17",
  "Th17", "ILC3", "CD8Tcm", "CD4Tcm", "M2", "NKreg", "PB", "Eosinophil",
  "PC", "cNK", "Tex", "intMo", "M1", "Tc"
)
heatmap_cells <- intersect(heatmap_cells, colnames(grouped_mean))
grouped_mean <- grouped_mean[, heatmap_cells, drop = FALSE]

message("Scaling each cell type across age groups...")
scaled_by_cell <- scale(grouped_mean)
scaled_by_cell[is.na(scaled_by_cell)] <- 0
scaled_by_cell[scaled_by_cell > 1] <- 1
scaled_by_cell[scaled_by_cell < -1] <- -1

age_category <- function(age_group) {
  ifelse(age_group %in% c("age0-1"), "infant_high",
         ifelse(age_group %in% c("age10-20", "age20-30", "age30-50"), "young_mid_high",
                "old_high"))
}

peak_index <- apply(scaled_by_cell, 2, which.max)
peak_group <- rownames(scaled_by_cell)[peak_index]
peak_value <- scaled_by_cell[cbind(peak_index, seq_along(peak_index))]
order_info <- data.frame(
  cell_type = colnames(scaled_by_cell),
  peak_group = peak_group,
  peak_category = age_category(peak_group),
  peak_value = as.numeric(peak_value),
  original_rank = match(colnames(scaled_by_cell), heatmap_cells),
  stringsAsFactors = FALSE
)
order_info$peak_category <- factor(order_info$peak_category,
                                   levels = c("infant_high", "young_mid_high", "old_high"))
order_info$peak_group <- factor(order_info$peak_group, levels = rownames(scaled_by_cell))
order_info <- order_info %>%
  arrange(peak_category, peak_group, desc(peak_value), original_rank)

display_names <- c(
  monoDC = "moDC",
  gdT = "γδT",
  Bnaive = "Bn"
)
ordered_cells <- order_info$cell_type
front_cells <- intersect(c("CD4Tn", "CD8Tn", "M0", "gdT"), ordered_cells)
ordered_cells <- c(front_cells, setdiff(ordered_cells, front_cells))
if (all(c("FOB", "CD8Temra") %in% ordered_cells)) {
  fob_i <- match("FOB", ordered_cells)
  cd8temra_i <- match("CD8Temra", ordered_cells)
  ordered_cells[c(fob_i, cd8temra_i)] <- ordered_cells[c(cd8temra_i, fob_i)]
}
order_info <- order_info[match(ordered_cells, order_info$cell_type), , drop = FALSE]
order_info$final_rank <- seq_len(nrow(order_info))
expr_set <- scaled_by_cell[, ordered_cells, drop = FALSE]
colnames(expr_set) <- ifelse(colnames(expr_set) %in% names(display_names),
                             display_names[colnames(expr_set)],
                             colnames(expr_set))
heatmap_mat <- expr_set

write.table(grouped_mean,
            file.path(out_dir, "ImmuCellAI2_age_group_mean_fraction_for_heatmap.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(heatmap_mat,
            file.path(out_dir, "ImmuCellAI2_age_group_zscore_heatmap_matrix.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(order_info,
            file.path(out_dir, "ImmuCellAI2_age_group_heatmap_cell_order.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Drawing heatmap...")
p_heatmap <- Heatmap(
  as.matrix(heatmap_mat),
  name = "Expression",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  width = ncol(heatmap_mat) * unit(4, "mm"),
  height = nrow(heatmap_mat) * unit(4, "mm"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = colorRampPalette(c("#3341A3", "black", "#E0DA54"))(50),
  row_dend_side = "left",
  column_dend_side = "top",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 8),
  border = TRUE,
  heatmap_legend_param = list(title = "Expression", at = c(-1, -0.5, 0, 0.5, 1),
                              labels = c("-1", "-0.5", "0", "0.5", "1"))
)

pdf_file <- file.path(out_dir, "Figure5_age_group_ImmuCellAI2_heatmap_Fig3E_like.pdf")
png_file <- file.path(out_dir, "Figure5_age_group_ImmuCellAI2_heatmap_Fig3E_like.png")

cairo_pdf(pdf_file, width = 12, height = 4.5)
draw(p_heatmap, heatmap_legend_side = "left", annotation_legend_side = "left",
     newpage = TRUE, merge_legend = TRUE)
dev.off()

png(png_file, width = 12, height = 4.5, units = "in", res = 600, type = "cairo")
draw(p_heatmap, heatmap_legend_side = "left", annotation_legend_side = "left",
     newpage = TRUE, merge_legend = TRUE)
dev.off()

cat("Done. Heatmap written to:\n", pdf_file, "\n", png_file, "\n", sep = "")
