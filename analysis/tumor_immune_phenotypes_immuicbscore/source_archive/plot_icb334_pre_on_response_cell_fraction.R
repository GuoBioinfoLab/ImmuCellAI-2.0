options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
})

base_dir <- "<LOCAL_R_ROOT>/Fig4"
out_dir <- file.path(base_dir, "ICB334_pre_on_response_cell_fraction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred_file <- file.path(
  base_dir,
  "Unified_ICB_score_cluster_ICB_integration",
  "fivefold_cv_unified_icb_score",
  "best_8cohort_pooled_5fold_predictions.txt"
)
frac_file <- file.path(
  base_dir,
  "C3_C4_score_external_immunotherapy_validation",
  "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt"
)
clin_file <- file.path(
  base_dir,
  "C3_C4_score_external_immunotherapy_validation",
  "external_immunotherapy_coldata_matched.txt"
)

target_cells <- c("cDC1", "Tc", "CD8Tem", "MBC")
cell_labels <- c(cDC1 = "cDC1", Tc = "Tc", CD8Tem = "CD8Tem", MBC = "Memory B cell")

theme_blue <- theme(
  plot.title = element_text(size = 13, face = "bold", color = "darkred", hjust = 0, lineheight = 1.2),
  plot.subtitle = element_text(size = 12, face = "bold", color = "grey30", lineheight = 1.2, hjust = 0),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = "gray80", linewidth = 0.7, linetype = "dashed"),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(vjust = 1, face = "bold", size = 12, color = "darkred"),
  axis.title.y = element_text(size = 12, face = "bold", color = "darkred"),
  axis.text.x = element_text(size = 11, colour = "black"),
  axis.text.y = element_text(size = 10, colour = "black"),
  legend.title = element_text(size = 11, colour = "black"),
  legend.text = element_text(size = 11, colour = "black"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
  legend.key = element_blank(),
  strip.background = element_rect(fill = "#F7E6E8", color = "black", linewidth = 0.6),
  strip.text = element_text(size = 11, colour = "black")
)

p_to_symbol <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

wilcox_one <- function(dat, group_col, value_col = "Cell_Ratio") {
  dat <- dat[!is.na(dat[[group_col]]) & !is.na(dat[[value_col]]), , drop = FALSE]
  groups_present <- unique(as.character(dat[[group_col]]))
  if (is.factor(dat[[group_col]])) {
    groups <- levels(dat[[group_col]])[levels(dat[[group_col]]) %in% groups_present]
  } else {
    groups <- sort(groups_present)
  }
  if (length(groups) != 2) {
    return(data.frame(group1 = NA_character_, group2 = NA_character_, p_value = NA_real_, statistic = NA_real_, n1 = NA_integer_, n2 = NA_integer_))
  }
  n1 <- sum(as.character(dat[[group_col]]) == groups[1])
  n2 <- sum(as.character(dat[[group_col]]) == groups[2])
  if (n1 < 2 || n2 < 2) {
    return(data.frame(group1 = groups[1], group2 = groups[2], p_value = NA_real_, statistic = NA_real_, n1 = n1, n2 = n2))
  }
  tmp <- dat
  tmp[[group_col]] <- factor(as.character(tmp[[group_col]]), levels = groups)
  wt <- suppressWarnings(wilcox.test(tmp[[value_col]] ~ tmp[[group_col]], exact = FALSE))
  data.frame(group1 = groups[1], group2 = groups[2], p_value = wt$p.value, statistic = unname(wt$statistic), n1 = n1, n2 = n2)
}

pred <- fread(pred_file)
frac <- fread(frac_file)
setnames(frac, 1, "sample")
clin <- fread(clin_file)

missing_cells <- setdiff(target_cells, colnames(frac))
if (length(missing_cells) > 0) {
  stop("Missing target cell columns: ", paste(missing_cells, collapse = ", "))
}

clin_keep <- clin %>%
  transmute(
    sample = Run,
    Clinical_Response = Response,
    Biopsy_Time_raw = Biopsy_Time,
    patient_id = patient_id,
    Response_standard = Response_standard
  )

plot_data <- pred %>%
  left_join(frac %>% select(sample, all_of(target_cells)), by = "sample") %>%
  left_join(clin_keep, by = "sample") %>%
  mutate(
    Response = case_when(
      ResponseBinary == 1 ~ "R",
      ResponseBinary == 0 ~ "NR",
      Clinical_Response %in% c("R", "NR") ~ Clinical_Response,
      TRUE ~ NA_character_
    ),
    Response = factor(Response, levels = c("NR", "R")),
    Biopsy_Time = case_when(
      grepl("^pre", Biopsy_Time_raw, ignore.case = TRUE) ~ "Pre",
      grepl("^on", Biopsy_Time_raw, ignore.case = TRUE) ~ "On",
      TRUE ~ NA_character_
    ),
    Biopsy_Time = factor(Biopsy_Time, levels = c("Pre", "On"))
  )

if (any(is.na(plot_data$Response)) || any(is.na(plot_data$Biopsy_Time))) {
  warning("Some samples have missing Response or Biopsy_Time after harmonization.")
}

long_data <- plot_data %>%
  select(sample, SRA_study, ResponseBinary, Unified_ICB_Score, Response, Biopsy_Time, all_of(target_cells)) %>%
  pivot_longer(cols = all_of(target_cells), names_to = "Cell_Type", values_to = "Cell_Ratio") %>%
  mutate(
    Cell_Type = factor(Cell_Type, levels = target_cells, labels = unname(cell_labels[target_cells])),
    Cell_Ratio = as.numeric(Cell_Ratio)
  )

fwrite(plot_data, file.path(out_dir, "ICB334_ImmuCellAI2_cell_fraction_with_clinical_wide.txt"), sep = "\t")
fwrite(long_data, file.path(out_dir, "ICB334_ImmuCellAI2_cell_fraction_with_clinical_long.txt"), sep = "\t")

sample_counts <- long_data %>%
  distinct(sample, SRA_study, Response, Biopsy_Time) %>%
  count(Response, Biopsy_Time, name = "N_samples")
fwrite(sample_counts, file.path(out_dir, "ICB334_response_time_sample_counts.txt"), sep = "\t")

time_stats <- long_data %>%
  group_by(Cell_Type, Response) %>%
  group_modify(~ wilcox_one(.x, "Biopsy_Time")) %>%
  ungroup() %>%
  mutate(comparison = "Pre_vs_On_within_Response", p_adj_BH = p.adjust(p_value, method = "BH"), signif = p_to_symbol(p_value))

response_on_stats <- long_data %>%
  filter(Biopsy_Time == "On") %>%
  group_by(Cell_Type) %>%
  group_modify(~ wilcox_one(.x, "Response")) %>%
  ungroup() %>%
  mutate(comparison = "R_vs_NR_within_On", p_adj_BH = p.adjust(p_value, method = "BH"), signif = p_to_symbol(p_value))

response_by_time_stats <- long_data %>%
  group_by(Cell_Type, Biopsy_Time) %>%
  group_modify(~ wilcox_one(.x, "Response")) %>%
  ungroup() %>%
  mutate(comparison = "R_vs_NR_within_Time", p_adj_BH = p.adjust(p_value, method = "BH"), signif = p_to_symbol(p_value))

all_stats <- bind_rows(
  time_stats,
  response_on_stats %>% mutate(Response = NA, Biopsy_Time = factor("On", levels = c("Pre", "On"))),
  response_by_time_stats %>% mutate(Response = NA)
)
fwrite(time_stats, file.path(out_dir, "ICB334_Pre_vs_On_within_response_wilcox.txt"), sep = "\t")
fwrite(response_on_stats, file.path(out_dir, "ICB334_On_treatment_R_vs_NR_wilcox.txt"), sep = "\t")
fwrite(response_by_time_stats, file.path(out_dir, "ICB334_R_vs_NR_within_each_time_wilcox.txt"), sep = "\t")
fwrite(all_stats, file.path(out_dir, "ICB334_all_requested_wilcox_statistics.txt"), sep = "\t")

response_cols <- c("NR" = "#4E79A7FF", "R" = "#E15759FF")
time_cols <- c("Pre" = "#76B7B2FF", "On" = "#F28E2BFF")

outlier_bounds <- long_data %>%
  group_by(Cell_Type, Response, Biopsy_Time) %>%
  summarise(
    q1 = as.numeric(quantile(Cell_Ratio, probs = 0.25, na.rm = TRUE)),
    q3 = as.numeric(quantile(Cell_Ratio, probs = 0.75, na.rm = TRUE)),
    iqr = IQR(Cell_Ratio, na.rm = TRUE),
    n_original = sum(!is.na(Cell_Ratio)),
    .groups = "drop"
  ) %>%
  mutate(
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr
  )

long_data_plot <- long_data %>%
  left_join(outlier_bounds, by = c("Cell_Type", "Response", "Biopsy_Time")) %>%
  mutate(is_plot_outlier = Cell_Ratio < lower | Cell_Ratio > upper) %>%
  filter(!is.na(Cell_Ratio), !is_plot_outlier) %>%
  mutate(Cell_Ratio_plot = Cell_Ratio)

outlier_counts <- long_data %>%
  left_join(outlier_bounds, by = c("Cell_Type", "Response", "Biopsy_Time")) %>%
  mutate(is_plot_outlier = Cell_Ratio < lower | Cell_Ratio > upper) %>%
  group_by(Cell_Type, Response, Biopsy_Time) %>%
  summarise(
    n_original = sum(!is.na(Cell_Ratio)),
    n_removed_for_plot = sum(is_plot_outlier, na.rm = TRUE),
    n_used_for_plot = n_original - n_removed_for_plot,
    lower = first(lower),
    upper = first(upper),
    .groups = "drop"
  )

plot_ranges_time <- long_data_plot %>%
  group_by(Cell_Type, Response) %>%
  summarise(
    y_label = max(Cell_Ratio_plot, na.rm = TRUE) * 0.94,
    .groups = "drop"
  ) %>%
  mutate(y_label = ifelse(is.finite(y_label) & y_label > 0, y_label, 1e-8))

plot_ranges_on <- long_data_plot %>%
  filter(Biopsy_Time == "On") %>%
  group_by(Cell_Type) %>%
  summarise(
    y_label = max(Cell_Ratio_plot, na.rm = TRUE) * 0.94,
    .groups = "drop"
  ) %>%
  mutate(y_label = ifelse(is.finite(y_label) & y_label > 0, y_label, 1e-8))

plot_ranges_response_by_time <- long_data_plot %>%
  group_by(Cell_Type, Biopsy_Time) %>%
  summarise(
    y_label = max(Cell_Ratio_plot, na.rm = TRUE) * 0.94,
    .groups = "drop"
  ) %>%
  mutate(y_label = ifelse(is.finite(y_label) & y_label > 0, y_label, 1e-8))

fwrite(outlier_counts, file.path(out_dir, "ICB334_plot_outlier_removed_counts.txt"), sep = "\t")

time_anno <- time_stats %>%
  left_join(plot_ranges_time, by = c("Cell_Type", "Response")) %>%
  mutate(x = 1.5, y = y_label)

response_on_anno <- response_on_stats %>%
  left_join(plot_ranges_on, by = "Cell_Type") %>%
  mutate(x = 1.5, y = y_label)

response_by_time_anno <- response_by_time_stats %>%
  left_join(plot_ranges_response_by_time, by = c("Cell_Type", "Biopsy_Time")) %>%
  mutate(x = 1.5, y = y_label)

p_time <- ggplot(long_data_plot, aes(Biopsy_Time, Cell_Ratio_plot, fill = Biopsy_Time)) +
  geom_violin(color = NA, alpha = 1, width = 0.72, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.42, linewidth = 0.55, alpha = 0.65, outlier.shape = NA) +
  facet_wrap(vars(Response, Cell_Type), scales = "free_y", ncol = 4) +
  scale_fill_manual(values = time_cols) +
  theme_blue +
  labs(
    title = "Pre- and on-treatment differences within response groups",
    subtitle = "Outliers are removed only for visualization by the 1.5IQR rule; Wilcoxon tests use original values",
    x = "",
    y = "Cell Ratio"
  ) +
  guides(fill = "none") +
  geom_text(
    data = time_anno,
    aes(x = x, y = y, label = signif),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  )

p_on_response <- ggplot(filter(long_data_plot, Biopsy_Time == "On"), aes(Response, Cell_Ratio_plot, fill = Response)) +
  geom_violin(color = NA, alpha = 1, width = 0.72, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.42, linewidth = 0.55, alpha = 0.65, outlier.shape = NA) +
  facet_wrap(~ Cell_Type, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = response_cols) +
  theme_blue +
  labs(
    title = "Response-group differences in on-treatment samples",
    subtitle = "Outliers are removed only for visualization by the 1.5IQR rule; Wilcoxon tests use original values",
    x = "",
    y = "Cell Ratio"
  ) +
  guides(fill = "none") +
  geom_text(
    data = response_on_anno,
    aes(x = x, y = y, label = signif),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  )

p_response_by_time <- ggplot(long_data_plot, aes(Response, Cell_Ratio_plot, fill = Response)) +
  geom_violin(color = NA, alpha = 1, width = 0.72, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.42, linewidth = 0.55, alpha = 0.65, outlier.shape = NA) +
  facet_wrap(vars(Biopsy_Time, Cell_Type), scales = "free_y", ncol = 4) +
  scale_fill_manual(values = response_cols) +
  theme_blue +
  labs(
    title = "Response-group differences stratified by biopsy time",
    subtitle = "Outliers are removed only for visualization by the 1.5IQR rule; Wilcoxon tests use original values",
    x = "",
    y = "Cell Ratio"
  ) +
  guides(fill = "none") +
  geom_text(
    data = response_by_time_anno,
    aes(x = x, y = y, label = signif),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  )

ggsave(file.path(out_dir, "Figure_ICB334_PreOn_within_Response_violin_boxplot.pdf"), p_time, width = 12, height = 6.8)
ggsave(file.path(out_dir, "Figure_ICB334_PreOn_within_Response_violin_boxplot.png"), p_time, width = 12, height = 6.8, dpi = 600)

ggsave(file.path(out_dir, "Figure_ICB334_OnTreatment_Response_violin_boxplot.pdf"), p_on_response, width = 12, height = 3.6)
ggsave(file.path(out_dir, "Figure_ICB334_OnTreatment_Response_violin_boxplot.png"), p_on_response, width = 12, height = 3.6, dpi = 600)

ggsave(file.path(out_dir, "Figure_ICB334_Response_by_BiopsyTime_violin_boxplot.pdf"), p_response_by_time, width = 12, height = 6.8)
ggsave(file.path(out_dir, "Figure_ICB334_Response_by_BiopsyTime_violin_boxplot.png"), p_response_by_time, width = 12, height = 6.8, dpi = 600)

p_combined <- plot_grid(
  p_time + theme(plot.title = element_text(size = 11, face = "bold", color = "darkred")),
  p_on_response + theme(plot.title = element_text(size = 11, face = "bold", color = "darkred")),
  ncol = 1,
  rel_heights = c(1, 0.42)
)

ggsave(file.path(out_dir, "Figure_ICB334_requested_comparison_combined.pdf"), p_combined, width = 12, height = 9)
ggsave(file.path(out_dir, "Figure_ICB334_requested_comparison_combined.png"), p_combined, width = 12, height = 9, dpi = 600)

cat("Output directory:", out_dir, "\n")
cat("N samples:", nrow(plot_data), "\n")
print(sample_counts)
cat("Done.\n")
