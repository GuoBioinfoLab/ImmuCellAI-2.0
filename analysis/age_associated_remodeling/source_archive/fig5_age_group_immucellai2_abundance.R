options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

fig5_dir <- "<LOCAL_R_ROOT>/Fig5"
out_dir <- file.path(fig5_dir, "ImmuCellAI2_age_group_abundance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bulk_file_1 <- file.path(fig5_dir, "age_bulk.txt")
bulk_file_2 <- file.path(fig5_dir, "age_bulk2.txt")
sample_info_file_1 <- file.path(fig5_dir, "age_sample_info.csv")
sample_info_file_2 <- file.path(fig5_dir, "healthy_info2_clean.csv")

reference_file <- "<LOCAL_R_ROOT>/reference_53celltypesTPM20260518.txt"
immune_gene_file <- "<LOCAL_R_ROOT>/MarkerUsedDeconvolution.txt"
package_file <- "<LOCAL_R_ROOT>/ImmuCellAI2.0_0.1.7.tar.gz"

stable_tmp <- file.path(out_dir, "Rtmp_immucellai2_age")
dir.create(stable_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = stable_tmp, TEMP = stable_tmp, TMPDIR = stable_tmp)

if (!requireNamespace("ImmuCellAI2.0", quietly = TRUE) && file.exists(package_file)) {
  install.packages(package_file, repos = NULL, type = "source")
}

suppressPackageStartupMessages(library(ImmuCellAI2.0))

n_cores <- min(8L, parallel::detectCores(logical = TRUE))

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

age_cols <- c(
  "age0-1" = "#E15759FF",
  "age10-20" = "#4E79A7FF",
  "age20-30" = "#F28E2BFF",
  "age30-50" = "#76B7B2FF",
  "age50-70" = "#59A14FFF",
  "age70+" = "#EDC948FF"
)

read_bulk_matrix <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  gene_col <- names(x)[1]
  genes <- as.character(x[[gene_col]])
  x[[gene_col]] <- NULL
  x[] <- lapply(x, function(v) as.numeric(as.character(v)))
  mat <- as.matrix(x)
  rownames(mat) <- genes
  mat <- mat[!is.na(rownames(mat)) & nzchar(rownames(mat)), , drop = FALSE]
  mat
}

collapse_duplicate_genes <- function(mat) {
  mat <- mat[!is.na(rownames(mat)) & nzchar(rownames(mat)), , drop = FALSE]
  if (!anyDuplicated(rownames(mat))) return(mat)
  collapsed <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  counts <- as.numeric(table(rownames(mat))[rownames(collapsed)])
  collapsed / counts
}

read_gene_list_local <- function(file) {
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- unlist(strsplit(x, "[,\t[:space:]]+"))
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!grepl("^[0-9]+$", x)]
  unique(x)
}

make_age_group <- function(age) {
  case_when(
    age >= 0 & age < 2 ~ "age0-1",
    age > 10 & age <= 20 ~ "age10-20",
    age > 20 & age <= 30 ~ "age20-30",
    age > 30 & age < 50 ~ "age30-50",
    age >= 50 & age < 70 ~ "age50-70",
    age >= 70 ~ "age70+",
    TRUE ~ NA_character_
  )
}

read_sample_info_1 <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  data.frame(
    Sample.ID = as.character(x[["Sample ID"]]),
    Age = as.numeric(x[["Age"]]),
    SourceBatch = "age_bulk",
    stringsAsFactors = FALSE
  )
}

read_sample_info_2 <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  data.frame(
    Sample.ID = as.character(x[["Run"]]),
    Age = as.numeric(x[["Age"]]),
    SourceBatch = "age_bulk2",
    RawGroup = if ("group" %in% names(x)) as.character(x[["group"]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

message("Reading bulk matrices...")
bulk1 <- collapse_duplicate_genes(read_bulk_matrix(bulk_file_1))
bulk2 <- collapse_duplicate_genes(read_bulk_matrix(bulk_file_2))

message("Reading sample age metadata...")
info1 <- read_sample_info_1(sample_info_file_1)
info2 <- read_sample_info_2(sample_info_file_2)
sample_info <- bind_rows(info1, info2) %>%
  mutate(
    Age = as.numeric(Age),
    group = make_age_group(Age),
    group = factor(group, levels = names(age_cols))
  ) %>%
  filter(!is.na(Sample.ID), !is.na(Age), !is.na(group))

bulk1_samples <- intersect(colnames(bulk1), sample_info$Sample.ID[sample_info$SourceBatch == "age_bulk"])
bulk2_samples <- intersect(colnames(bulk2), sample_info$Sample.ID[sample_info$SourceBatch == "age_bulk2"])
bulk1 <- bulk1[, bulk1_samples, drop = FALSE]
bulk2 <- bulk2[, bulk2_samples, drop = FALSE]

common_bulk_genes <- intersect(rownames(bulk1), rownames(bulk2))
if (length(common_bulk_genes) < 1000L) {
  stop("Too few common genes between age_bulk.txt and age_bulk2.txt: ", length(common_bulk_genes))
}
bulk <- cbind(bulk1[common_bulk_genes, , drop = FALSE], bulk2[common_bulk_genes, , drop = FALSE])

sample_info <- sample_info %>%
  filter(Sample.ID %in% colnames(bulk)) %>%
  distinct(Sample.ID, .keep_all = TRUE)
sample_order <- intersect(colnames(bulk), sample_info$Sample.ID)
bulk <- bulk[, sample_order, drop = FALSE]
sample_info <- sample_info[match(sample_order, sample_info$Sample.ID), , drop = FALSE]

message("Reading reference matrix...")
reference <- collapse_duplicate_genes(read_expression_matrix(reference_file))
common_genes <- intersect(rownames(bulk), rownames(reference))
bulk <- bulk[common_genes, , drop = FALSE]
reference <- reference[common_genes, , drop = FALSE]

keep <- rowSums(bulk, na.rm = TRUE) > 0 & rowSums(reference, na.rm = TRUE) > 0
bulk <- bulk[keep, , drop = FALSE]
reference <- reference[keep, , drop = FALSE]

immune_gene_tokens <- read_gene_list_local(immune_gene_file)
immune_genes <- intersect(immune_gene_tokens, rownames(reference))
if (length(immune_genes) >= 100L) {
  bulk <- bulk[immune_genes, , drop = FALSE]
  reference <- reference[immune_genes, , drop = FALSE]
}

write.table(sample_info, file.path(out_dir, "ImmuCellAI2_age_combined_sample_info_with_groups.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

group_counts <- sample_info %>%
  count(group, SourceBatch, name = "N_samples") %>%
  arrange(group, SourceBatch)
write.table(group_counts, file.path(out_dir, "ImmuCellAI2_age_group_sample_counts_by_batch.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

qc <- data.frame(
  Item = c("age_bulk_samples", "age_bulk2_samples", "combined_samples", "common_bulk_genes",
           "genes_after_reference_intersection", "marker_gene_tokens",
           "marker_genes_overlap_reference_bulk", "genes_used_after_marker_filter"),
  Value = c(length(bulk1_samples), length(bulk2_samples), ncol(bulk), length(common_bulk_genes),
            length(common_genes), length(immune_gene_tokens), length(immune_genes), nrow(bulk))
)
write.table(qc, file.path(out_dir, "ImmuCellAI2_age_input_QC.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

state_fraction_file <- file.path(out_dir, "ImmuCellAI2_age_state_fraction_sample_by_celltype.txt")
if (file.exists(state_fraction_file)) {
  message("Using existing ImmuCellAI2 state fraction file: ", state_fraction_file)
  state_fraction <- read.table(state_fraction_file, header = TRUE, sep = "\t", row.names = 1,
                               check.names = FALSE, comment.char = "")
  state_fraction <- as.matrix(state_fraction)
  state_fraction <- state_fraction[intersect(sample_order, rownames(state_fraction)), , drop = FALSE]
  sample_order <- rownames(state_fraction)
  sample_info <- sample_info[match(sample_order, sample_info$Sample.ID), , drop = FALSE]
} else {
  message("Running ImmuCellAI2 VB deconvolution...")
  hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = "tcell")
  hierarchy <- ImmuCellAI2.0:::validate_hierarchy(hierarchy, reference)
  reference <- reference[, hierarchy$state_name, drop = FALSE]

  fit <- deconvolve_bulk_matrix(
    bulk.mat = bulk,
    reference = reference,
    hierarchy = hierarchy,
    hierarchy.mode = "tcell",
    inference.method = "vb",
    add.unknown = FALSE,
    pseudo.depth = 1e5,
    n.iter = 50,
    vb.tol = 1e-6,
    alpha.major = 10,
    alpha.sub = 5,
    alpha.state = 1,
    n.cores = n_cores,
    seed = 123
  )

  state_fraction <- as.matrix(fit$state.fraction)
  state_fraction <- state_fraction[sample_order, , drop = FALSE]

  write.table(state_fraction, state_fraction_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(t(state_fraction), file.path(out_dir, "ImmuCellAI2_age_state_fraction_celltype_by_sample.txt"),
              sep = "\t", quote = FALSE, col.names = NA)
}

plot_data <- as.data.frame(state_fraction) %>%
  tibble::rownames_to_column("Sample.ID") %>%
  left_join(sample_info %>% select(Sample.ID, Age, group, SourceBatch), by = "Sample.ID") %>%
  pivot_longer(cols = all_of(colnames(state_fraction)), names_to = "cell_type", values_to = "value") %>%
  mutate(
    cell_type = recode(cell_type, Bnaive = "Bn", monoDC = "moDC"),
    value = as.numeric(value)
  )

write.table(plot_data, file.path(out_dir, "ImmuCellAI2_age_long_data_for_plot.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

flag_plot_outliers <- function(dat) {
  dat %>%
    group_by(cell_type, group) %>%
    mutate(
      plot_q1 = quantile(value, 0.25, na.rm = TRUE),
      plot_q3 = quantile(value, 0.75, na.rm = TRUE),
      plot_iqr = plot_q3 - plot_q1,
      plot_lower = plot_q1 - 1.5 * plot_iqr,
      plot_upper = plot_q3 + 1.5 * plot_iqr,
      plot_outlier = if_else(
        is.na(value) | is.na(plot_iqr) | plot_iqr <= 0,
        FALSE,
        value < plot_lower | value > plot_upper
      )
    ) %>%
    ungroup()
}

plot_data_with_outlier_flag <- flag_plot_outliers(plot_data)
plot_data_no_outliers <- plot_data_with_outlier_flag %>%
  filter(!plot_outlier | cell_type == "CD8Tem") %>%
  select(-plot_q1, -plot_q3, -plot_iqr, -plot_lower, -plot_upper)

write.table(plot_data_with_outlier_flag, file.path(out_dir, "ImmuCellAI2_age_long_data_for_plot_with_outlier_flag.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(plot_data_no_outliers, file.path(out_dir, "ImmuCellAI2_age_long_data_for_plot_no_outliers.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

group_counts_total <- plot_data %>%
  distinct(Sample.ID, Age, group, SourceBatch) %>%
  count(group, name = "N_samples")
write.table(group_counts_total, file.path(out_dir, "ImmuCellAI2_age_group_sample_counts.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cell_y_limits <- list(
  M0 = c(0, 0.018),
  Bn = c(0, 0.01),
  CD8Tem = c(0, 0.0038)
)
outlier_keep_cells <- c("CD4Tem", "CD8Tem")
outlier_force_remove_cells <- c("CD8Temra")

make_cell_plot <- function(cell, jitter_size = 1.8, remove_outliers = TRUE) {
  use_filtered_data <- (remove_outliers || cell %in% outlier_force_remove_cells) &&
    !(cell %in% outlier_keep_cells)
  dat <- if (use_filtered_data) {
    plot_data_no_outliers %>% filter(cell_type == cell)
  } else {
    plot_data %>% filter(cell_type == cell)
  }
  p <- ggplot(dat, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.color = NA, width = 0.8) +
    geom_jitter(shape = 21, color = "gray2", position = position_jitter(0.17), size = jitter_size, alpha = 0.78) +
    scale_x_discrete(labels = c("0-1", "10-20", "20-30", "30-50", "50-70", "70+")) +
    scale_fill_manual(values = age_cols, drop = FALSE) +
    labs(title = cell, x = "Age Group", y = "Cell Ratio", fill = NULL) +
    theme_blue +
    geom_smooth(
      method = "loess", formula = y ~ x, linewidth = 1, se = FALSE,
      linetype = "solid", aes(group = 1), color = "#3F88C5"
    ) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  if (cell %in% names(cell_y_limits)) {
    p <- p + coord_cartesian(ylim = cell_y_limits[[cell]])
  }
  p
}

save_cell_panel <- function(cells, figure_name, ncol, width, height, jitter_size = 1.8, dpi = 600, remove_outliers = TRUE) {
  missing_cells <- setdiff(cells, unique(plot_data$cell_type))
  if (length(missing_cells) > 0) {
    warning("Missing cells skipped in ", figure_name, ": ", paste(missing_cells, collapse = ", "))
  }
  cells <- intersect(cells, unique(plot_data$cell_type))
  panel_plots <- lapply(cells, make_cell_plot, jitter_size = jitter_size, remove_outliers = remove_outliers)
  names(panel_plots) <- cells
  p <- plot_grid(plotlist = panel_plots, ncol = ncol, align = "hv")
  ggsave(file.path(out_dir, paste0(figure_name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(out_dir, paste0(figure_name, ".png")), p, width = width, height = height, dpi = dpi)
  invisible(p)
}

fig3a_cells <- c("CD4Tn", "CD8Tn", "M0", "gdT")
fig3b_cells <- c("cNK", "CD8Temra", "ncMo", "Eosinophil")
fig3c_cells <- c("CD4Tem", "CD8Tem", "CD4Tcm", "Basophil", "cMo", "Breg", "Bn", "MZB", "ILC1")

p_fig3a <- save_cell_panel(fig3a_cells, "Figure5_age_group_ImmuCellAI2_Fig3A_like_cells", ncol = 4, width = 16, height = 4.4, remove_outliers = FALSE)
p_fig3b <- save_cell_panel(fig3b_cells, "Figure5_age_group_ImmuCellAI2_Fig3B_like_cells", ncol = 4, width = 16, height = 4.4, remove_outliers = FALSE)
p_fig3c <- save_cell_panel(fig3c_cells, "Figure5_age_group_ImmuCellAI2_Fig3C_like_cells", ncol = 3, width = 12, height = 13.2, remove_outliers = TRUE)

p_fig3c_combined <- plot_grid(p_fig3c, NULL, ncol = 2, rel_widths = c(3, 1))
p_selected <- plot_grid(p_fig3a, p_fig3b, p_fig3c_combined, ncol = 1, rel_heights = c(1, 1, 3))
ggsave(file.path(out_dir, "Figure5_age_group_ImmuCellAI2_Fig3ABC_like_combined.pdf"), p_selected, width = 18, height = 22)
ggsave(file.path(out_dir, "Figure5_age_group_ImmuCellAI2_Fig3ABC_like_combined.png"), p_selected, width = 18, height = 22, dpi = 500)

all53_cells <- colnames(state_fraction)
all53_cells <- recode(all53_cells, Bnaive = "Bn", monoDC = "moDC")
all53_cells <- all53_cells[all53_cells %in% unique(plot_data$cell_type)]
p_all53 <- save_cell_panel(
  all53_cells,
  "Figure5_age_group_ImmuCellAI2_all_53_cells",
  ncol = 6,
  width = 24,
  height = 34,
  jitter_size = 1.15,
  dpi = 450
)

trend_summary <- plot_data %>%
  group_by(cell_type) %>%
  summarise(
    N = sum(!is.na(value) & !is.na(Age)),
    spearman_age = suppressWarnings(cor(Age, value, method = "spearman", use = "complete.obs")),
    pearson_age = suppressWarnings(cor(Age, value, method = "pearson", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(spearman_age)))
write.table(trend_summary, file.path(out_dir, "ImmuCellAI2_age_celltype_correlation_with_age.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

run_info <- data.frame(
  Item = c("InputBulk1", "InputBulk2", "SampleInfo1", "SampleInfo2", "Reference", "ImmuneGeneFile",
           "Samples", "GenesUsed", "HierarchyMode", "InferenceMethod", "AddUnknown", "Cores", "PackageVersion"),
  Value = c(bulk_file_1, bulk_file_2, sample_info_file_1, sample_info_file_2, reference_file, immune_gene_file,
            ncol(bulk), nrow(bulk), "tcell", "vb", "FALSE", n_cores,
            as.character(packageVersion("ImmuCellAI2.0")))
)
write.table(run_info, file.path(out_dir, "ImmuCellAI2_age_run_info.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done. Results written to:", out_dir, "\n")
print(group_counts_total)
