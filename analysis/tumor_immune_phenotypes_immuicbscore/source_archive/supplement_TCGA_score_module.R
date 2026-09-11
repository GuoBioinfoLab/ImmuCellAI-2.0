options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(cowplot)
})

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
source_dir <- file.path(fig4_dir, "cluster4_class_Unified_ICB_score_association_full")
extended_dir <- file.path(fig4_dir, "Unified_ICB_score_cluster_ICB_integration", "extended_TCGA_cluster_score_analysis")
out_dir <- file.path(fig4_dir, "TCGA_score_supplement_module")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(source_dir, "cluster4_class_samples_with_unified_ICB_score.txt")
quartile_file <- file.path(source_dir, "cluster4_class_composition_by_score_quartile.txt")
or_file <- file.path(extended_dir, "logistic_OR_per_1SD_score_for_cluster_membership.txt")
effect_file <- file.path(extended_dir, "pairwise_cluster_score_effect_size.txt")
decile_file <- file.path(extended_dir, "score_decile_cluster_composition.txt")

stopifnot(file.exists(score_file), file.exists(quartile_file), file.exists(or_file), file.exists(effect_file), file.exists(decile_file))

cluster_cols <- c(
  "C1" = "#E15759FF",
  "C2" = "#4E79A7FF",
  "C3" = "#F28E2BFF",
  "C4" = "#76B7B2FF"
)

theme_supp <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey86", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      strip.background = element_rect(fill = "#F7E7EA", color = "black", linewidth = 0.45),
      strip.text = element_text(color = "black", face = "bold"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold", color = "darkred", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(color = "grey30", hjust = 0, size = base_size)
    )
}

save_plot <- function(p, name, width, height) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height, units = "in")
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, units = "in", dpi = 500)
}

score_df <- fread(score_file, data.table = FALSE) %>%
  mutate(
    cluster = factor(cluster, levels = paste0("C", 1:4)),
    score_quartile = factor(score_quartile, levels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"))
  )

quartile_df <- fread(quartile_file, data.table = FALSE) %>%
  mutate(
    cluster = factor(cluster, levels = paste0("C", 1:4)),
    score_quartile = factor(score_quartile, levels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"))
  )

or_df <- fread(or_file, data.table = FALSE) %>%
  mutate(cluster = factor(cluster, levels = paste0("C", 1:4)))

effect_df <- fread(effect_file, data.table = FALSE) %>%
  filter(cluster_b == "C4" | cluster_a == "C4") %>%
  mutate(
    comparison_to_C4 = ifelse(cluster_a == "C4", paste0("C4 vs ", cluster_b), paste0("C4 vs ", cluster_a)),
    cliffs_delta_C4_minus_other = ifelse(cluster_a == "C4", cliffs_delta_a_minus_b, -cliffs_delta_a_minus_b),
    mean_diff_C4_minus_other = ifelse(cluster_a == "C4", mean_diff_a_minus_b, -mean_diff_a_minus_b),
    comparison_to_C4 = factor(comparison_to_C4, levels = c("C4 vs C1", "C4 vs C2", "C4 vs C3"))
  )

decile_df <- fread(decile_file, data.table = FALSE) %>%
  mutate(
    cluster = factor(cluster, levels = paste0("C", 1:4)),
    score_decile = factor(score_decile, levels = paste0("D", 1:10))
  )

summary_tbl <- score_df %>%
  group_by(cluster) %>%
  summarise(
    N = n(),
    mean_score = mean(Unified_ICB_score, na.rm = TRUE),
    median_score = median(Unified_ICB_score, na.rm = TRUE),
    sd_score = sd(Unified_ICB_score, na.rm = TRUE),
    .groups = "drop"
  )
fwrite(summary_tbl, file.path(out_dir, "TCGA_cluster_ImmuICBscore_summary.txt"), sep = "\t")
fwrite(or_df, file.path(out_dir, "TCGA_cluster_membership_OR_per_1SD_ImmuICBscore.txt"), sep = "\t")
fwrite(effect_df, file.path(out_dir, "TCGA_C4_vs_other_cluster_ImmuICBscore_effect_size.txt"), sep = "\t")

p_a <- ggplot(score_df, aes(x = cluster, y = Unified_ICB_score, fill = cluster)) +
  geom_violin(width = 0.86, color = "black", linewidth = 0.3, trim = TRUE) +
  geom_boxplot(width = 0.16, color = "black", linewidth = 0.3, outlier.shape = NA, alpha = 0.82) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
  scale_fill_manual(values = cluster_cols) +
  labs(
    title = "A. ImmuICBscore in TCGA immune clusters",
    x = "TCGA immune cluster",
    y = "ImmuICBscore"
  ) +
  theme_supp() +
  theme(legend.position = "none")

p_b <- ggplot(quartile_df, aes(x = score_quartile, y = Proportion, fill = cluster)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "B. Cluster composition across score quartiles",
    x = "ImmuICBscore quartile",
    y = "Cluster proportion",
    fill = "Cluster"
  ) +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_c <- ggplot(or_df, aes(x = OR_per_1SD_score, y = cluster, color = cluster)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_errorbar(aes(xmin = Lower95, xmax = Upper95), orientation = "y", width = 0.18, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_color_manual(values = cluster_cols) +
  scale_x_log10() +
  labs(
    title = "C. Cluster membership odds per 1SD score",
    x = "Odds ratio per 1SD higher ImmuICBscore",
    y = "TCGA immune cluster"
  ) +
  theme_supp() +
  theme(legend.position = "none")

p_d <- ggplot(effect_df, aes(x = comparison_to_C4, y = cliffs_delta_C4_minus_other, fill = comparison_to_C4)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.35) +
  geom_col(width = 0.62, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("C4 vs C1" = "#76B7B2FF", "C4 vs C2" = "#59A14FFF", "C4 vs C3" = "#B07AA1FF")) +
  labs(
    title = "D. C4 score advantage over other clusters",
    x = "",
    y = "Cliff's delta of ImmuICBscore"
  ) +
  theme_supp() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

p_e <- ggplot(decile_df, aes(x = score_decile, y = Proportion, group = cluster, color = cluster)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(values = cluster_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "E. Cluster trends across score deciles",
    x = "ImmuICBscore decile",
    y = "Cluster proportion",
    color = "Cluster"
  ) +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_a, "Supplement_TCGA_score_module_A_score_by_cluster", 4.5, 3.6)
save_plot(p_b, "Supplement_TCGA_score_module_B_quartile_composition", 5.2, 3.6)
save_plot(p_c, "Supplement_TCGA_score_module_C_OR_per_1SD", 4.5, 3.4)
save_plot(p_d, "Supplement_TCGA_score_module_D_C4_effect_size", 4.3, 3.4)
save_plot(p_e, "Supplement_TCGA_score_module_E_decile_trend", 6.2, 3.6)

top_row <- plot_grid(p_a, p_b, ncol = 2, rel_widths = c(0.9, 1.1), labels = NULL)
mid_row <- plot_grid(p_c, p_d, ncol = 2, rel_widths = c(1, 1), labels = NULL)
p_combined <- plot_grid(top_row, mid_row, p_e, ncol = 1, rel_heights = c(1, 1, 1.05))

save_plot(p_combined, "Supplement_TCGA_ImmuICBscore_relationship_module_combined", 9.2, 10.5)

pdf(file.path(out_dir, "Supplement_TCGA_ImmuICBscore_relationship_module_multipage.pdf"), width = 7.2, height = 5)
print(p_a)
print(p_b)
print(p_c)
print(p_d)
print(p_e)
dev.off()

readme <- c(
  "Supplementary module: TCGA immune clusters and ImmuICBscore relationship",
  "",
  "Purpose:",
  "This module shows that the immune-cell-composition-based ICB response score is associated with TCGA immune cluster structure.",
  "",
  "Panels:",
  "A. ImmuICBscore distribution across TCGA C1-C4 immune clusters.",
  "B. TCGA cluster composition across ImmuICBscore quartiles.",
  "C. Logistic odds ratio for membership in each cluster per 1SD increase of ImmuICBscore.",
  "D. C4 versus other clusters: score effect size using Cliff's delta.",
  "E. Cluster proportion trends across ImmuICBscore deciles.",
  "",
  "Interpretation:",
  "C4 shows the highest mean ImmuICBscore and the strongest positive OR per 1SD score. C1 is depleted as score increases, supporting the relationship between the score and the TCGA immune-cluster axis.",
  "",
  "Input files:",
  score_file,
  quartile_file,
  or_file,
  effect_file,
  decile_file
)
writeLines(readme, file.path(out_dir, "README_TCGA_ImmuICBscore_supplement_module.txt"))

cat("Supplementary TCGA-score module written to:", out_dir, "\n")
