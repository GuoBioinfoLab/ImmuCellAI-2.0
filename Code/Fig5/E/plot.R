# E is a data-derived heatmap, not a manually drawn illustration.
heatmap_mat <- analysis$E
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
plot_files <- c(plot_files, save_plot(p_heatmap, "E", width = 12, height = 4.5, heatmap = TRUE))
