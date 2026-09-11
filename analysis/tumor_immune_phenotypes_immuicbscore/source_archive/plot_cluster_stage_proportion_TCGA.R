options(stringsAsFactors = FALSE)

cluster_file <- "<LOCAL_CLUSTER_ROOT>/cluster_TCGA.k=4.consensusClass.csv"
clinical_file <- "<LOCAL_CLUSTER_ROOT>/clinical_TCGA.txt"
immucellai2_file <- "<LOCAL_R_ROOT>/Fig4/TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000/TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt"
out_dir <- "<LOCAL_R_ROOT>/Fig4"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

normalize_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("^'--$|^--$|^NA$|^NOT REPORTED$|^NOT APPLICABLE$|^UNKNOWN$|^$", NA_character_, x)
  x <- gsub("STAGE\\s+", "STAGE ", x)
  out <- rep(NA_character_, length(x))
  out[grepl("^STAGE\\s+I($|[^I,V])|^I($|[^I,V])", x)] <- "Stage I"
  out[grepl("^STAGE\\s+II($|[^I])|^II($|[^I])", x)] <- "Stage II"
  out[grepl("^STAGE\\s+III|^III", x)] <- "Stage III"
  out[grepl("^STAGE\\s+IV|^IV", x)] <- "Stage IV"
  factor(out, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
}

message("Reading cluster class...")
cluster4_class <- fread(cluster_file, header = FALSE, data.table = FALSE)
colnames(cluster4_class) <- c("sample", "cluster")
cluster4_class$sample <- as.character(cluster4_class$sample)
cluster4_class$case_submitter_id <- substr(cluster4_class$sample, 1, 12)
cluster4_class$cluster <- factor(
  paste0("Cluster", as.integer(cluster4_class$cluster)),
  levels = paste0("Cluster", 1:4)
)

message("Reading ImmuCellAI2 result sample names...")
immune_samples <- fread(immucellai2_file, nrows = 1, data.table = FALSE)
result_TCGA <- fread(immucellai2_file, select = 1, data.table = FALSE)
immune_sample_names <- as.character(result_TCGA[[1]])

cluster4_class <- cluster4_class %>%
  filter(sample %in% immune_sample_names)

message("Reading clinical stage...")
clinical_lines <- readLines(clinical_file, warn = FALSE, encoding = "UTF-8")
clinical_lines <- clinical_lines[-1]
clinical_match <- regexec(
  "^\\S+\\s+\\S+\\s+(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})\\s+(TCGA-[A-Z0-9]+)\\s+\\S+\\s+((?:Stage\\s+[IVX]+[A-Z]*)|'--|--)",
  clinical_lines,
  ignore.case = FALSE
)
clinical_parts <- regmatches(clinical_lines, clinical_match)
clinical_ok <- lengths(clinical_parts) == 4L
if (!any(clinical_ok)) {
  stop("Could not parse clinical stage rows from: ", clinical_file)
}
clinical <- data.frame(
  case_submitter_id = vapply(clinical_parts[clinical_ok], `[`, character(1), 2),
  project_id = vapply(clinical_parts[clinical_ok], `[`, character(1), 3),
  ajcc_pathologic_stage = vapply(clinical_parts[clinical_ok], `[`, character(1), 4),
  stringsAsFactors = FALSE
)

clinical_stage_long <- clinical %>%
  transmute(
    case_submitter_id = as.character(case_submitter_id),
    project_id = as.character(project_id),
    ajcc_pathologic_stage_raw = as.character(ajcc_pathologic_stage),
    Stage = normalize_stage(ajcc_pathologic_stage)
  ) %>%
  filter(!is.na(Stage))

stage_conflict <- clinical_stage_long %>%
  distinct(case_submitter_id, Stage) %>%
  count(case_submitter_id, name = "N_stage") %>%
  filter(N_stage > 1)

write.table(
  stage_conflict,
  file.path(out_dir, "TCGA_cluster4_stage_conflicting_cases.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# For duplicate clinical records, keep the most frequent normalized stage per case.
clinical_stage <- clinical_stage_long %>%
  count(case_submitter_id, project_id, Stage, name = "N_records") %>%
  group_by(case_submitter_id) %>%
  arrange(desc(N_records), Stage, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

merged <- cluster4_class %>%
  left_join(clinical_stage, by = "case_submitter_id") %>%
  filter(!is.na(Stage))

stage_counts <- merged %>%
  count(cluster, Stage, name = "N") %>%
  group_by(cluster) %>%
  mutate(
    ClusterN = sum(N),
    Proportion = N / ClusterN
  ) %>%
  ungroup() %>%
  tidyr::complete(
    cluster = factor(paste0("Cluster", 1:4), levels = paste0("Cluster", 1:4)),
    Stage = factor(c("Stage I", "Stage II", "Stage III", "Stage IV"),
                   levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    fill = list(N = 0, Proportion = 0)
  ) %>%
  group_by(cluster) %>%
  mutate(ClusterN = sum(N)) %>%
  ungroup()

cluster_totals <- merged %>%
  count(cluster, name = "N_with_stage")

write.table(
  merged,
  file.path(out_dir, "TCGA_cluster4_samples_with_stage.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  stage_counts,
  file.path(out_dir, "TCGA_cluster4_stage_proportion_table.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  cluster_totals,
  file.path(out_dir, "TCGA_cluster4_stage_sample_counts.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

stage_cols <- c(
  "Stage I" = "#E15759FF",
  "Stage II" = "#4E79A7FF",
  "Stage III" = "#F28E2BFF",
  "Stage IV" = "#76B7B2FF"
)

p <- ggplot(stage_counts, aes(x = cluster, y = Proportion, fill = Stage)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(Proportion >= 0.04, percent(Proportion, accuracy = 1), "")),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  scale_fill_manual(values = stage_cols, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Clinical stage proportion", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45)
  )

pdf(file.path(out_dir, "TCGA_cluster4_stage_proportion_stacked_bar.pdf"),
    width = 6.2, height = 4.8, useDingbats = FALSE)
print(p)
dev.off()

png(file.path(out_dir, "TCGA_cluster4_stage_proportion_stacked_bar.png"),
    width = 2200, height = 1600, res = 300)
print(p)
dev.off()

message("Samples in cluster result and ImmuCellAI2 result: ", nrow(cluster4_class))
message("Samples with Stage I-IV information: ", nrow(merged))
message("Conflicting clinical stage cases: ", nrow(stage_conflict))
print(stage_counts)
message("Done. Outputs written to: ", out_dir)
