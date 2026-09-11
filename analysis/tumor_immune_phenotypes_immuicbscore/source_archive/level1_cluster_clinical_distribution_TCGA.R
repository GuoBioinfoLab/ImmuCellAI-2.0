options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_file <- "<LOCAL_CLUSTER_ROOT>/cluster_TCGA.k=4.consensusClass.csv"
sample_info_file <- file.path(fig4_dir, "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000", "TCGA_ImmunotherapyResponsive_tumor_sample_info.txt")
immucellai2_file <- file.path(fig4_dir, "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000", "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
clinical_file <- file.path(fig4_dir, "clinical.tsv")
out_dir <- file.path(fig4_dir, "level1_cluster_clinical_distribution")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

clean_value <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "'--", "--", "NA", "N/A", "na", "not reported", "not applicable", "unknown", "Unknown")] <- NA_character_
  x
}

mode_value <- function(x) {
  x <- clean_value(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

median_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(clean_value(x)))
  if (all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

normalize_stage <- function(pathologic_stage, clinical_stage = NA_character_) {
  x <- clean_value(pathologic_stage)
  y <- clean_value(clinical_stage)
  x[is.na(x)] <- y[is.na(x)]
  xu <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(xu))
  out[grepl("^STAGE\\s+I($|[^I,V])|^I($|[^I,V])", xu)] <- "Stage I"
  out[grepl("^STAGE\\s+II($|[^I])|^II($|[^I])", xu)] <- "Stage II"
  out[grepl("^STAGE\\s+III|^III", xu)] <- "Stage III"
  out[grepl("^STAGE\\s+IV|^IV", xu)] <- "Stage IV"
  factor(out, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
}

normalize_grade <- function(x) {
  x <- toupper(clean_value(x))
  out <- rep(NA_character_, length(x))
  out[grepl("G1|GRADE 1|LOW GRADE", x)] <- "Grade 1"
  out[grepl("G2|GRADE 2", x)] <- "Grade 2"
  out[grepl("G3|GRADE 3|HIGH GRADE", x)] <- "Grade 3"
  out[grepl("G4|GRADE 4", x)] <- "Grade 4"
  out[is.na(out) & !is.na(x)] <- x[is.na(out) & !is.na(x)]
  factor(out, levels = c("Grade 1", "Grade 2", "Grade 3", "Grade 4", sort(setdiff(unique(out), c("Grade 1", "Grade 2", "Grade 3", "Grade 4", NA)))))
}

categorical_test <- function(dat, variable) {
  sub <- dat[, c("cluster", variable)]
  sub <- sub[complete.cases(sub), , drop = FALSE]
  if (nrow(sub) == 0 || length(unique(sub[[variable]])) < 2 || length(unique(sub$cluster)) < 2) {
    return(data.frame(Variable = variable, Test = NA_character_, P_value = NA_real_, N = nrow(sub), stringsAsFactors = FALSE))
  }
  tab <- table(sub$cluster, sub[[variable]])
  p <- tryCatch({
    suppressWarnings(chisq.test(tab)$p.value)
  }, error = function(e) NA_real_)
  data.frame(Variable = variable, Test = "Chi-square", P_value = p, N = nrow(sub), stringsAsFactors = FALSE)
}

continuous_test <- function(dat, variable) {
  sub <- dat[, c("cluster", variable)]
  sub <- sub[complete.cases(sub), , drop = FALSE]
  if (nrow(sub) == 0 || length(unique(sub$cluster)) < 2) {
    return(data.frame(Variable = variable, Test = NA_character_, P_value = NA_real_, N = nrow(sub), stringsAsFactors = FALSE))
  }
  p <- tryCatch({
    kruskal.test(sub[[variable]] ~ sub$cluster)$p.value
  }, error = function(e) NA_real_)
  data.frame(Variable = variable, Test = "Kruskal-Wallis", P_value = p, N = nrow(sub), stringsAsFactors = FALSE)
}

make_prop_table <- function(dat, variable) {
  sub <- dat[, c("cluster", variable)]
  colnames(sub) <- c("cluster", "Feature")
  sub <- sub[complete.cases(sub), , drop = FALSE]
  if (!nrow(sub)) {
    return(data.frame(cluster = character(), Feature = character(), N = integer(), ClusterN = integer(), Proportion = numeric()))
  }
  sub %>%
    count(cluster, Feature, name = "N") %>%
    group_by(cluster) %>%
    mutate(ClusterN = sum(N), Proportion = N / ClusterN) %>%
    ungroup()
}

plot_stacked_prop <- function(prop_df, variable, fill_values = NULL, width = 7.2, height = 4.8) {
  if (!nrow(prop_df)) return(invisible(NULL))
  p <- ggplot(prop_df, aes(x = cluster, y = Proportion, fill = Feature)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.25) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Proportion", fill = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 11, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45)
    )
  if (!is.null(fill_values)) {
    p <- p + scale_fill_manual(values = fill_values, drop = FALSE)
  }
  pdf(file.path(out_dir, paste0("TCGA_cluster4_", variable, "_stacked_bar.pdf")), width = width, height = height, useDingbats = FALSE)
  print(p)
  dev.off()
  png(file.path(out_dir, paste0("TCGA_cluster4_", variable, "_stacked_bar.png")), width = width * 300, height = height * 300, res = 300)
  print(p)
  dev.off()
  invisible(p)
}

message("Reading cluster, sample info, and ImmuCellAI2 sample names...")
cluster4 <- fread(cluster_file, header = FALSE, data.table = FALSE)
colnames(cluster4) <- c("Sample", "ClusterNumber")
cluster4$Sample <- as.character(cluster4$Sample)
cluster4$CaseSubmitterID <- substr(cluster4$Sample, 1, 12)
cluster4$cluster <- factor(paste0("Cluster", as.integer(cluster4$ClusterNumber)), levels = paste0("Cluster", 1:4))

sample_info <- fread(sample_info_file, data.table = FALSE)
immune_samples <- fread(immucellai2_file, select = 1, data.table = FALSE)[[1]]

cluster4 <- cluster4 %>%
  filter(Sample %in% immune_samples) %>%
  left_join(sample_info, by = c("Sample", "CaseSubmitterID"))

message("Reading clinical.tsv...")
clinical <- fread(clinical_file, sep = "\t", data.table = FALSE)

clinical_case <- clinical %>%
  group_by(case_submitter_id) %>%
  summarise(
    project_id_clinical = mode_value(project_id),
    age_at_index = median_numeric(age_at_index),
    age_at_diagnosis_days = median_numeric(age_at_diagnosis),
    gender = mode_value(gender),
    race = mode_value(race),
    ethnicity = mode_value(ethnicity),
    vital_status = mode_value(vital_status),
    ajcc_pathologic_stage_raw = mode_value(ajcc_pathologic_stage),
    ajcc_clinical_stage_raw = mode_value(ajcc_clinical_stage),
    tumor_grade_raw = mode_value(tumor_grade),
    prior_malignancy = mode_value(prior_malignancy),
    prior_treatment = mode_value(prior_treatment),
    progression_or_recurrence = mode_value(progression_or_recurrence),
    .groups = "drop"
  ) %>%
  mutate(
    Stage = normalize_stage(ajcc_pathologic_stage_raw, ajcc_clinical_stage_raw),
    TumorGrade = normalize_grade(tumor_grade_raw),
    age_at_diagnosis_years = age_at_diagnosis_days / 365.25
  )

analysis_df <- cluster4 %>%
  left_join(clinical_case, by = c("CaseSubmitterID" = "case_submitter_id")) %>%
  mutate(
    ProjectID = ifelse(is.na(ProjectID), project_id_clinical, ProjectID),
    ProjectID = factor(ProjectID),
    SampleType = factor(SampleType),
    gender = factor(gender),
    race = factor(race),
    ethnicity = factor(ethnicity),
    vital_status = factor(vital_status),
    prior_malignancy = factor(prior_malignancy),
    prior_treatment = factor(prior_treatment),
    progression_or_recurrence = factor(progression_or_recurrence)
  )

write.table(analysis_df, file.path(out_dir, "TCGA_cluster4_level1_clinical_merged_samples.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

categorical_vars <- c(
  "ProjectID", "Stage", "gender", "SampleType", "race", "ethnicity",
  "vital_status", "TumorGrade", "prior_malignancy", "prior_treatment",
  "progression_or_recurrence"
)

prop_tables <- list()
for (v in categorical_vars) {
  prop <- make_prop_table(analysis_df, v)
  if (nrow(prop) > 0) {
    prop_tables[[v]] <- data.frame(Variable = v, prop, stringsAsFactors = FALSE)
  } else {
    prop_tables[[v]] <- data.frame(
      Variable = character(),
      cluster = character(),
      Feature = character(),
      N = integer(),
      ClusterN = integer(),
      Proportion = numeric(),
      stringsAsFactors = FALSE
    )
  }
  write.table(prop, file.path(out_dir, paste0("TCGA_cluster4_", v, "_proportion_table.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
all_prop <- bind_rows(prop_tables)
write.table(all_prop, file.path(out_dir, "TCGA_cluster4_all_categorical_proportion_tables.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

test_results <- bind_rows(
  lapply(categorical_vars, function(v) categorical_test(analysis_df, v)),
  continuous_test(analysis_df, "age_at_index"),
  continuous_test(analysis_df, "age_at_diagnosis_years")
) %>%
  mutate(FDR = p.adjust(P_value, method = "BH"))
write.table(test_results, file.path(out_dir, "TCGA_cluster4_level1_clinical_association_tests.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cluster_summary <- analysis_df %>%
  group_by(cluster) %>%
  summarise(
    N = n(),
    N_stage = sum(!is.na(Stage)),
    N_age = sum(!is.na(age_at_index)),
    MedianAge = median(age_at_index, na.rm = TRUE),
    MeanAge = mean(age_at_index, na.rm = TRUE),
    .groups = "drop"
  )
write.table(cluster_summary, file.path(out_dir, "TCGA_cluster4_level1_cluster_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

stage_cols <- c("Stage I" = "#E15759FF", "Stage II" = "#4E79A7FF", "Stage III" = "#F28E2BFF", "Stage IV" = "#76B7B2FF")
gender_cols <- c("female" = "#E15759FF", "male" = "#4E79A7FF")
sample_type_cols <- c("Primary Tumor" = "#4E79A7FF", "Recurrent Solid Tumor" = "#F28E2BFF", "Additional New Primary" = "#59A14FFF", "Metastatic" = "#E15759FF")

project_palette <- setNames(
  grDevices::colorRampPalette(c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF", "#59A14FFF", "#EDC948FF", "#B07AA1FF", "#FF9DA7FF", "#9C755FFF", "#BAB0ACFF"))(length(unique(na.omit(analysis_df$ProjectID)))),
  sort(unique(as.character(na.omit(analysis_df$ProjectID))))
)

plot_stacked_prop(prop_tables$ProjectID, "ProjectID", project_palette, width = 8.5, height = 5.2)
plot_stacked_prop(prop_tables$Stage, "Stage", stage_cols, width = 6.2, height = 4.8)
plot_stacked_prop(prop_tables$gender, "gender", gender_cols, width = 5.6, height = 4.5)
plot_stacked_prop(prop_tables$SampleType, "SampleType", sample_type_cols, width = 6.2, height = 4.8)
plot_stacked_prop(prop_tables$TumorGrade, "TumorGrade", NULL, width = 6.2, height = 4.8)
plot_stacked_prop(prop_tables$vital_status, "vital_status", NULL, width = 5.8, height = 4.6)

age_df <- analysis_df %>% filter(!is.na(age_at_index))
p_age <- ggplot(age_df, aes(x = cluster, y = age_at_index, fill = cluster)) +
  geom_boxplot(width = 0.62, outlier.size = 0.4, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("#E15759FF", "#4E79A7FF", "#F28E2BFF", "#76B7B2FF"), guide = "none") +
  labs(x = NULL, y = "Age at index") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45)
  )
pdf(file.path(out_dir, "TCGA_cluster4_age_at_index_boxplot.pdf"), width = 5.6, height = 4.5, useDingbats = FALSE)
print(p_age)
dev.off()
png(file.path(out_dir, "TCGA_cluster4_age_at_index_boxplot.png"), width = 1680, height = 1350, res = 300)
print(p_age)
dev.off()

message("Done. Outputs written to: ", out_dir)
print(cluster_summary)
print(test_results)
