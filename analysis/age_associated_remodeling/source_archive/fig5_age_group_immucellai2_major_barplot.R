options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

fig5_dir <- "<LOCAL_R_ROOT>/Fig5"
out_dir <- file.path(fig5_dir, "ImmuCellAI2_age_group_abundance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fraction_file <- file.path(out_dir, "ImmuCellAI2_age_state_fraction_sample_by_celltype.txt")
sample_info_file <- file.path(out_dir, "ImmuCellAI2_age_combined_sample_info_with_groups.txt")

age_levels <- c("age0-1", "age10-20", "age20-30", "age30-50", "age50-70", "age70+")

theme_blue <- theme(
  plot.title = element_text(size = 13, face = "bold", color = "darkred", hjust = 0, lineheight = 1.2),
  plot.subtitle = element_text(size = 11, face = "bold", color = "grey30", lineheight = 1.2, hjust = 0),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = "gray80", linewidth = 0.7, linetype = "dashed"),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(vjust = 1, face = "bold", size = 12, color = "darkred"),
  axis.title.y = element_text(size = 12, face = "bold", color = "darkred"),
  axis.text.x = element_text(size = 10, colour = "black"),
  axis.text.y = element_text(size = 10, colour = "black"),
  legend.title = element_text(size = 11, colour = "black"),
  legend.text = element_text(size = 10, colour = "black"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
  legend.key = element_blank(),
  strip.background = element_rect(fill = "#F7E6E8", color = "black", linewidth = 0.6),
  strip.text = element_text(size = 10, colour = "black", face = "bold")
)

read_fraction <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  sample_col <- names(x)[1]
  rownames(x) <- x[[sample_col]]
  x[[sample_col]] <- NULL
  x[] <- lapply(x, function(v) as.numeric(as.character(v)))
  as.data.frame(x, check.names = FALSE)
}

sum_present <- function(dat, cells) {
  cells <- intersect(cells, colnames(dat))
  if (length(cells) == 0) return(rep(0, nrow(dat)))
  rowSums(dat[, cells, drop = FALSE], na.rm = TRUE)
}

message("Reading ImmuCellAI2 age fractions...")
fraction <- read_fraction(fraction_file)

message("Reading age metadata...")
sample_info <- fread(sample_info_file, data.table = FALSE, check.names = FALSE) %>%
  mutate(group = factor(group, levels = age_levels)) %>%
  filter(Sample.ID %in% rownames(fraction), !is.na(group))

fraction <- fraction[sample_info$Sample.ID, , drop = FALSE]

major_fraction <- data.frame(
  Sample.ID = sample_info$Sample.ID,
  group = sample_info$group,
  Bcells = sum_present(fraction, c("Bnaive", "Breg", "FOB", "BGC", "MZB", "MBC", "PC", "PB", "Bex")),
  Tcells = sum_present(fraction, c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm",
                                   "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm",
                                   "MAIT", "NKT", "Tc", "Tex", "Tfh", "Th1", "Th1/17",
                                   "Th17", "Th2", "Tr1", "Treg", "gdT")),
  ILCs = sum_present(fraction, c("ILC1", "ILC2", "ILC3")),
  NK = sum_present(fraction, c("NKreg", "cNK")),
  Macrophages = sum_present(fraction, c("M0", "M1", "M2", "TAM")),
  Granulocytes = sum_present(fraction, c("Neutrophil", "Basophil", "Eosinophil", "Mast cell")),
  DCs = sum_present(fraction, c("pDC", "cDC1", "cDC2", "monoDC", "Langerhans")),
  Monocytes = sum_present(fraction, c("cMo", "intMo", "ncMo", "MDSC")),
  stringsAsFactors = FALSE
)

major_fraction$major_sum <- rowSums(major_fraction[, c("Bcells", "Tcells", "ILCs", "NK",
                                                       "Macrophages", "Granulocytes", "DCs",
                                                       "Monocytes")], na.rm = TRUE)

major_group_mean <- major_fraction %>%
  group_by(group) %>%
  summarise(across(c(Bcells, Tcells, ILCs, NK, Macrophages, Granulocytes, DCs, Monocytes, major_sum),
                   function(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(group = factor(group, levels = age_levels))

plot_data <- major_group_mean %>%
  select(-major_sum) %>%
  pivot_longer(cols = -group, names_to = "cell_type", values_to = "proportion")

cell_order <- c("Tcells", "Bcells", "ILCs", "NK",
                "Granulocytes", "Monocytes", "Macrophages", "DCs")
plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_order)
plot_data$group <- factor(plot_data$group, levels = age_levels)

cell_colors <- c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF",
                 "#59A14FFF", "#EDC948FF", "#B07AA1FF", "#FF9DA7FF")
names(cell_colors) <- cell_order

write.table(major_fraction,
            file.path(out_dir, "ImmuCellAI2_age_major_cell_fraction_sample_level.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(major_group_mean,
            file.path(out_dir, "ImmuCellAI2_age_major_cell_fraction_group_mean.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(plot_data,
            file.path(out_dir, "ImmuCellAI2_age_major_cell_fraction_barplot_long.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Drawing major cell class stacked bar plot...")
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

pdf_file <- file.path(out_dir, "Figure5_age_group_ImmuCellAI2_major_cell_barplot_Fig3D_like.pdf")
png_file <- file.path(out_dir, "Figure5_age_group_ImmuCellAI2_major_cell_barplot_Fig3D_like.png")

ggsave(pdf_file, pD, width = 5, height = 7)
ggsave(png_file, pD, width = 5, height = 7, dpi = 600)

cat("Done. Major cell bar plot written to:\n", pdf_file, "\n", png_file, "\n", sep = "")
