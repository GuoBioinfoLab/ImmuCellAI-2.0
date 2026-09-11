options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

fig6_dir <- "<LOCAL_R_ROOT>/Fig6"
out_dir <- file.path(fig6_dir, "ImmuCellAI2_active_latent_TB")
metadata_file <- file.path(fig6_dir, "PRJNA_combined.csv")
state_fraction_file <- file.path(out_dir, "ImmuCellAI2_TB_state_fraction_sample_by_celltype.txt")

metadata <- fread(metadata_file, data.table = FALSE, check.names = FALSE) %>%
  filter(tb_status %in% c("Active TB", "Latent TB")) %>%
  distinct(Run, .keep_all = TRUE)
metadata$tb_status <- factor(metadata$tb_status, levels = c("Active TB", "Latent TB"))

state_fraction <- read.table(
  state_fraction_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  comment.char = ""
)

selected_map <- c(Tex = "Tex", PB = "PB")
missing_cells <- setdiff(unname(selected_map), colnames(state_fraction))
if (length(missing_cells) > 0) {
  stop("Missing cells in ImmuCellAI2 result: ", paste(missing_cells, collapse = ", "))
}

common_samples <- intersect(metadata$Run, rownames(state_fraction))
metadata <- metadata[match(common_samples, metadata$Run), , drop = FALSE]
state_fraction <- as.matrix(state_fraction[common_samples, unname(selected_map), drop = FALSE])

plot_df <- as.data.frame(state_fraction) %>%
  tibble::rownames_to_column("Run")
colnames(plot_df)[match(unname(selected_map), colnames(plot_df))] <- names(selected_map)

plot_df <- plot_df %>%
  left_join(metadata %>% select(Run, BioProject, Disease, tb_status), by = "Run") %>%
  pivot_longer(cols = all_of(names(selected_map)), names_to = "CellType", values_to = "Abundance") %>%
  mutate(
    status = factor(tb_status, levels = c("Active TB", "Latent TB")),
    CellType = factor(CellType, levels = c("Tex", "PB")),
    Abundance = as.numeric(Abundance)
  )

stats <- plot_df %>%
  group_by(CellType) %>%
  summarise(
    N_Active = sum(status == "Active TB" & !is.na(Abundance)),
    N_Latent = sum(status == "Latent TB" & !is.na(Abundance)),
    Mean_Active = mean(Abundance[status == "Active TB"], na.rm = TRUE),
    Mean_Latent = mean(Abundance[status == "Latent TB"], na.rm = TRUE),
    Median_Active = median(Abundance[status == "Active TB"], na.rm = TRUE),
    Median_Latent = median(Abundance[status == "Latent TB"], na.rm = TRUE),
    FoldChange_Latent_vs_Active = Mean_Latent / Mean_Active,
    P_Wilcox = suppressWarnings(wilcox.test(Abundance ~ status)$p.value),
    .groups = "drop"
  ) %>%
  mutate(
    Direction = case_when(
      Mean_Latent > Mean_Active ~ "Latent TB higher",
      Mean_Latent < Mean_Active ~ "Active TB higher",
      TRUE ~ "No mean difference"
    )
  )

write.table(
  plot_df,
  file.path(out_dir, "ImmuCellAI2_TB_Tex_PB_Active_vs_Latent_long.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  stats,
  file.path(out_dir, "ImmuCellAI2_TB_Tex_PB_Active_vs_Latent_stats.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sig_symbol <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

make_half_violin_data <- function(dat, width = 0.31) {
  dat <- dat %>% filter(!is.na(Abundance), !is.na(status), !is.na(CellType))
  status_levels <- levels(dat$status)
  pieces <- list()
  k <- 1L
  for (status_i in status_levels) {
    vals <- dat$Abundance[dat$status == status_i]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L || length(unique(vals)) < 2L) next
    dens <- density(vals, n = 128, from = min(vals), to = max(vals), na.rm = TRUE)
    dens_scaled <- dens$y / max(dens$y, na.rm = TRUE) * width
    x_base <- ifelse(status_i == "Active TB", 1.07, 2.07)
    pieces[[k]] <- data.frame(
      status = factor(status_i, levels = status_levels),
      x = c(x_base + dens_scaled, rep(x_base, length(dens$x))),
      y = c(dens$x, rev(dens$x)),
      group_id = status_i,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  bind_rows(pieces)
}

make_panel <- function(cell, show_y = FALSE) {
  tmp <- plot_df %>%
    filter(CellType == cell) %>%
    mutate(
      x = ifelse(status == "Active TB", 1, 2),
      x_box = x - 0.08
    )
  half_df <- make_half_violin_data(tmp)
  pval <- stats$P_Wilcox[as.character(stats$CellType) == cell]
  y_max <- max(tmp$Abundance, na.rm = TRUE)
  y_min <- min(tmp$Abundance, na.rm = TRUE)
  y_pad <- ifelse(y_max > y_min, (y_max - y_min) * 0.18, max(y_max, 1e-6) * 0.25)
  y_label <- y_max + y_pad * 0.45

  ggplot(tmp, aes(x = x_box, y = Abundance)) +
    geom_polygon(
      data = half_df,
      aes(x = x, y = y, group = group_id, fill = status),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.8,
      alpha = 0.9
    ) +
    stat_boxplot(
      aes(group = status, colour = status),
      geom = "errorbar",
      width = 0.18,
      linewidth = 0.75
    ) +
    geom_boxplot(
      aes(group = status, colour = status, fill = status),
      width = 0.34,
      outlier.shape = NA,
      linewidth = 0.75,
      alpha = 0.95
    ) +
    stat_summary(
      aes(group = status),
      fun = "median",
      geom = "point",
      shape = 16,
      size = 2.1,
      color = "white"
    ) +
    annotate("text", x = 1.5, y = y_label, label = sig_symbol(pval),
             size = 8, fontface = "bold", color = "black") +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Active TB", "Latent TB"),
      limits = c(0.55, 2.55),
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = c("Active TB" = "#E15759FF", "Latent TB" = "#4E79A7FF")) +
    scale_color_manual(values = c("Active TB" = "black", "Latent TB" = "black")) +
    labs(x = NULL, y = if (show_y) "Cell Ratio" else NULL, title = cell) +
    coord_cartesian(ylim = c(y_min, y_max + y_pad), clip = "off") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 15, face = "bold", color = "#8B1A1A", hjust = 0),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.y = element_text(size = 12, color = "#8B1A1A", face = "bold"),
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.8, linetype = "dashed"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.05),
      aspect.ratio = 1,
      plot.margin = margin(4, 8, 4, 8)
    )
}

p <- cowplot::plot_grid(
  make_panel("Tex", TRUE),
  make_panel("PB", FALSE),
  nrow = 1,
  align = "hv",
  labels = c("A", ""),
  label_size = 18,
  label_fontface = "bold",
  hjust = -0.2,
  vjust = 1.4
)

pdf_file <- file.path(out_dir, "Figure6_Tex_PB_ImmuCellAI2_ActiveTB_vs_LatentTB.pdf")
png_file <- file.path(out_dir, "Figure6_Tex_PB_ImmuCellAI2_ActiveTB_vs_LatentTB.png")
ggsave(pdf_file, p, width = 6.4, height = 3.2, bg = "white")
ggsave(png_file, p, width = 6.4, height = 3.2, dpi = 600, bg = "white")

cat("Done. Tex/PB comparison written to:\n", out_dir, "\n", sep = "")
print(stats)
