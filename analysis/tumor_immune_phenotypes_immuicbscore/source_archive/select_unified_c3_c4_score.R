options(stringsAsFactors = FALSE)

base_dir <- "<LOCAL_R_ROOT>/Fig4/C3_C4_score_external_immunotherapy_validation/C3_C4_axis_ICB_response_scores"
score_file <- file.path(base_dir, "external_ICB_C3_C4_axis_scores.txt")
meta_file <- "<LOCAL_R_ROOT>/Fig4/C3_C4_score_external_immunotherapy_validation/external_immunotherapy_C3_C4_scores_with_clinical.txt"

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

auc_manual <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

scores <- fread(score_file, data.table = FALSE, check.names = FALSE)
meta <- fread(meta_file, data.table = FALSE, check.names = FALSE)

response <- rep(NA_integer_, nrow(meta))
response[meta$ResponseGroup == "Responder"] <- 1L
response[meta$ResponseGroup == "Non-responder"] <- 0L
response[is.na(response) & meta$Response == "R"] <- 1L
response[is.na(response) & meta$Response == "NR"] <- 0L
meta$response_binary <- response
meta <- meta[!is.na(meta$response_binary) & meta$Run %in% scores$sample, , drop = FALSE]
scores <- scores[match(meta$Run, scores$sample), , drop = FALSE]

candidate_scores <- data.frame(
  sample = scores$sample,
  C3_low_clr = -scores$C3_pearson_clr,
  C3_low_spearman = -scores$C3_spearman_raw,
  C4vC3_ranger_asin = scores$TCGA_C4_vs_C3_ranger_probability_asin,
  C4vC3_ranger_clr = scores$TCGA_C4_vs_C3_ranger_probability_clr,
  C4minusC3_cosine = scores$C4_minus_C3_cosine_raw,
  Top5_C4_over_C3 = scores$auto_top5_C4sum_over_C3sum_logratio,
  stringsAsFactors = FALSE
)

rows <- list()
for (score_name in setdiff(colnames(candidate_scores), "sample")) {
  score_value <- candidate_scores[[score_name]]
  rows[[length(rows) + 1L]] <- data.frame(
    score = score_name,
    subgroup_variable = "ALL",
    subgroup = "ALL",
    N = length(score_value),
    N_R = sum(meta$response_binary == 1),
    N_NR = sum(meta$response_binary == 0),
    AUC = auc_manual(meta$response_binary, score_value),
    stringsAsFactors = FALSE
  )

  for (subgroup_variable in c("SRA_study", "Cancer", "Cancer_type", "Drug", "Anti_target", "disease")) {
    if (!subgroup_variable %in% colnames(meta)) next
    vals <- as.character(meta[[subgroup_variable]])
    vals[is.na(vals) | vals == ""] <- "Unknown"
    for (subgroup in sort(unique(vals))) {
      idx <- which(vals == subgroup)
      y <- meta$response_binary[idx]
      if (length(idx) < 20 || sum(y == 1) < 5 || sum(y == 0) < 5) next
      rows[[length(rows) + 1L]] <- data.frame(
        score = score_name,
        subgroup_variable = subgroup_variable,
        subgroup = subgroup,
        N = length(idx),
        N_R = sum(y == 1),
        N_NR = sum(y == 0),
        AUC = auc_manual(y, score_value[idx]),
        stringsAsFactors = FALSE
      )
    }
  }
}

auc_by_subgroup <- rbindlist(rows)

sra_summary <- auc_by_subgroup[subgroup_variable == "SRA_study", .(
  MeanAUC = mean(AUC, na.rm = TRUE),
  MedianAUC = median(AUC, na.rm = TRUE),
  WeightedAUC = weighted.mean(AUC, N, na.rm = TRUE),
  MinAUC = min(AUC, na.rm = TRUE),
  MaxAUC = max(AUC, na.rm = TRUE),
  N_subgroups = .N,
  N_AUC_ge_0.7 = sum(AUC >= 0.7, na.rm = TRUE),
  N_AUC_ge_0.8 = sum(AUC >= 0.8, na.rm = TRUE)
), by = score][order(-MeanAUC)]

all_summary <- auc_by_subgroup[subgroup_variable == "ALL"][order(-AUC)]

fwrite(candidate_scores, file.path(base_dir, "unified_candidate_scores_for_ICB_samples.txt"), sep = "\t")
fwrite(auc_by_subgroup, file.path(base_dir, "unified_score_candidate_auc_by_subgroup_fixed_direction.txt"), sep = "\t")
fwrite(sra_summary, file.path(base_dir, "unified_score_candidate_auc_summary_by_SRA_study_fixed_direction.txt"), sep = "\t")
fwrite(all_summary, file.path(base_dir, "unified_score_candidate_auc_overall_fixed_direction.txt"), sep = "\t")

final_score <- "C3_low_clr"
final_df <- data.frame(
  Run = meta$Run,
  Response = meta$Response,
  ResponseGroup = meta$ResponseGroup,
  response_binary = meta$response_binary,
  SRA_study = if ("SRA_study" %in% colnames(meta)) meta$SRA_study else NA_character_,
  Cancer = if ("Cancer" %in% colnames(meta)) meta$Cancer else NA_character_,
  Cancer_type = if ("Cancer_type" %in% colnames(meta)) meta$Cancer_type else NA_character_,
  Drug = if ("Drug" %in% colnames(meta)) meta$Drug else NA_character_,
  disease = if ("disease" %in% colnames(meta)) meta$disease else NA_character_,
  Unified_C3low_ICB_score = candidate_scores[[final_score]],
  stringsAsFactors = FALSE
)
fwrite(final_df, file.path(base_dir, "final_unified_C3low_ICB_score_with_clinical.txt"), sep = "\t")

plot_df <- auc_by_subgroup[subgroup_variable == "SRA_study" & score == final_score]
plot_df <- plot_df[order(-plot_df$AUC), ]
p <- ggplot(plot_df, aes(x = reorder(subgroup, AUC), y = AUC)) +
  geom_col(fill = "#379DA5", color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0.8, color = "#F66463", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "grey40", linetype = "dotted") +
  coord_flip() +
  theme_bw(base_size = 9) +
  labs(
    x = NULL,
    y = "AUC",
    title = "Unified C3-low ICB score across external SRA-study cohorts",
    subtitle = "Score = -cor(CLR(sample cell fractions), CLR(TCGA C3 centroid))"
  ) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(base_dir, "final_unified_C3low_ICB_score_SRA_study_auc.pdf"), p, width = 7.2, height = 5)
ggsave(file.path(base_dir, "final_unified_C3low_ICB_score_SRA_study_auc.png"), p, width = 7.2, height = 5, dpi = 300)

print(sra_summary)
print(all_summary)
print(plot_df)
