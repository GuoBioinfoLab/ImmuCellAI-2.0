plot_files <- c(plot_files, local({
  summary_long <- prepare_plot_data("E")$plot_data
  summary_long$Time <- factor(summary_long$Time, levels = c("D0", "D1", "D4", "D8"))
  phy_cols <- c("Died" = "#A3C8DC", "Survived" = "#EA5D2D")
p <- ggplot(summary_long, aes(x = Time, y = Mean, group = status, color = status)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ CellType, ncol = 4, scales = "free_y") +
  labs(x = "Time", y = "Ratios", title = "Cell Type Ratios by Status") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.75, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.05),
    aspect.ratio = 1,
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 14, color = "black", hjust = 0),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  scale_colour_manual(values = phy_cols, drop = FALSE)


  save_plot(p, "E")
}))
