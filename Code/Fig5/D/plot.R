# Sourced by plot_panels.R; no A-C display filtering is applied to D.
cell_order <- c("Tcells", "Bcells", "ILCs", "NK",
                "Granulocytes", "Monocytes", "Macrophages", "DCs")
cell_colors <- c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF",
                 "#59A14FFF", "#EDC948FF", "#B07AA1FF", "#FF9DA7FF")
names(cell_colors) <- cell_order
plot_data <- analysis$D
plot_data$group <- factor(plot_data$group, levels = age_levels)
plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_order)
pD <- ggplot(plot_data, aes(x = group, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.78, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = cell_colors, drop = FALSE) +
  labs(x = "Age Group", y = "Frequency (%)", title = "Cell Distribution by Age") +
  theme_blue +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 100, 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "on") +
  scale_x_discrete(labels = c("0-1", "10-20", "20-30", "30-50", "50-70", "70+"), drop = FALSE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    plot.title = element_text(size = 13, face = "bold", color = "darkred")
  )
plot_files <- c(plot_files, save_plot(pD, "D", width = 5, height = 7))
