options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

fig6_dir <- "<LOCAL_R_ROOT>/Fig6"
out_dir <- file.path(fig6_dir, "ImmuCellAI2_PRJNA683803_AIDS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

exp_file <- file.path(fig6_dir, "exp.txt")
metadata_file <- file.path(fig6_dir, "PRJNA683803_combined.csv")
reference_file <- "<LOCAL_R_ROOT>/reference_53celltypesTPM20260518.txt"
immune_gene_file <- "<LOCAL_R_ROOT>/MarkerUsedDeconvolution.txt"
package_file <- "<LOCAL_R_ROOT>/ImmuCellAI2.0_0.1.7.tar.gz"

stable_tmp <- file.path(out_dir, "Rtmp_immucellai2_aids")
dir.create(stable_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = stable_tmp, TEMP = stable_tmp, TMPDIR = stable_tmp)

if (!requireNamespace("ImmuCellAI2.0", quietly = TRUE) && file.exists(package_file)) {
  install.packages(package_file, repos = NULL, type = "source")
}
suppressPackageStartupMessages(library(ImmuCellAI2.0))

n_cores <- min(8L, parallel::detectCores(logical = TRUE))

read_gene_list_local <- function(file) {
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- unlist(strsplit(x, "[,\t[:space:]]+"))
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!grepl("^[0-9]+$", x)]
  unique(x)
}

collapse_duplicate_genes <- function(mat) {
  mat <- mat[!is.na(rownames(mat)) & nzchar(rownames(mat)), , drop = FALSE]
  if (!anyDuplicated(rownames(mat))) return(mat)
  collapsed <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  counts <- as.numeric(table(rownames(mat))[rownames(collapsed)])
  collapsed / counts
}

normalize_outcome <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    tolower(x) %in% c("died", "dead") ~ "Died",
    tolower(x) %in% c("survived", "survive") ~ "Survived",
    tolower(x) %in% c("died-iris", "dead-iris") ~ "Died-IRIS",
    tolower(x) %in% c("survived-iris", "survive-iris") ~ "Survived-IRIS",
    TRUE ~ x
  )
}

normalize_treatment <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    tolower(x) %in% c("early", "early art", "earlyart") ~ "Early ART",
    tolower(x) %in% c("deferred", "deferred art", "deferredart") ~ "Deferred ART",
    TRUE ~ x
  )
}

message("Reading PRJNA683803 metadata...")
metadata <- fread(metadata_file, data.table = FALSE, check.names = FALSE) %>%
  distinct(Run, .keep_all = TRUE) %>%
  mutate(
    outcome = normalize_outcome(outcome),
    treatment = normalize_treatment(treatment)
  )

message("Checking expression columns...")
exp_header <- names(fread(exp_file, nrows = 0, data.table = FALSE, check.names = FALSE))
matched_samples <- metadata$Run[metadata$Run %in% exp_header]
if (length(matched_samples) == 0L) {
  stop("No PRJNA683803 metadata samples are present in exp.txt.")
}
metadata <- metadata[match(matched_samples, metadata$Run), , drop = FALSE]

state_fraction_file <- file.path(out_dir, "ImmuCellAI2_PRJNA683803_state_fraction_sample_by_celltype.txt")
if (file.exists(state_fraction_file)) {
  message("Using existing ImmuCellAI2 result: ", state_fraction_file)
  state_fraction <- read.table(
    state_fraction_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    comment.char = ""
  )
  state_fraction <- as.matrix(state_fraction)
  state_fraction <- state_fraction[matched_samples, , drop = FALSE]
  genes_used <- NA_integer_
  marker_tokens_n <- NA_integer_
  marker_overlap_n <- NA_integer_
} else {
  message("Reading selected expression columns: ", length(matched_samples), " samples...")
  exp_df <- fread(
    exp_file,
    data.table = FALSE,
    check.names = FALSE,
    select = c("V1", matched_samples)
  )
  gene_col <- names(exp_df)[1]
  genes <- as.character(exp_df[[gene_col]])
  exp_df[[gene_col]] <- NULL
  exp_df[] <- lapply(exp_df, function(v) as.numeric(as.character(v)))
  bulk <- as.matrix(exp_df)
  rownames(bulk) <- genes
  bulk <- collapse_duplicate_genes(bulk)
  bulk <- bulk[, matched_samples, drop = FALSE]

  message("Reading reference and filtering marker genes...")
  reference <- collapse_duplicate_genes(read_expression_matrix(reference_file))
  common_genes <- intersect(rownames(bulk), rownames(reference))
  bulk <- bulk[common_genes, , drop = FALSE]
  reference <- reference[common_genes, , drop = FALSE]

  keep <- rowSums(bulk, na.rm = TRUE) > 0 & rowSums(reference, na.rm = TRUE) > 0
  bulk <- bulk[keep, , drop = FALSE]
  reference <- reference[keep, , drop = FALSE]

  marker_tokens <- read_gene_list_local(immune_gene_file)
  marker_genes <- intersect(marker_tokens, rownames(reference))
  if (length(marker_genes) >= 100L) {
    bulk <- bulk[marker_genes, , drop = FALSE]
    reference <- reference[marker_genes, , drop = FALSE]
  }

  hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = "tcell")
  hierarchy <- ImmuCellAI2.0:::validate_hierarchy(hierarchy, reference)
  reference <- reference[, hierarchy$state_name, drop = FALSE]

  message("Running ImmuCellAI2 VB deconvolution for PRJNA683803...")
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
  state_fraction <- state_fraction[matched_samples, , drop = FALSE]
  write.table(state_fraction, state_fraction_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(
    t(state_fraction),
    file.path(out_dir, "ImmuCellAI2_PRJNA683803_state_fraction_celltype_by_sample.txt"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  genes_used <- nrow(bulk)
  marker_tokens_n <- length(marker_tokens)
  marker_overlap_n <- length(marker_genes)
}

write.table(
  metadata,
  file.path(out_dir, "PRJNA683803_matched_metadata.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

state_df <- as.data.frame(state_fraction) %>%
  tibble::rownames_to_column("Run") %>%
  left_join(metadata %>% select(Run, BioProject, Disease, outcome, treatment), by = "Run")

write.table(
  state_df,
  file.path(out_dir, "ImmuCellAI2_PRJNA683803_state_fraction_with_metadata.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mark_iqr_outlier <- function(x) {
  ok <- is.finite(x)
  out <- rep(FALSE, length(x))
  if (sum(ok) < 4L || length(unique(x[ok])) < 2L) return(out)
  qs <- quantile(x[ok], probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  iqr <- qs[2] - qs[1]
  if (!is.finite(iqr) || iqr <= 0) return(out)
  out[ok] <- x[ok] < (qs[1] - 1.5 * iqr) | x[ok] > (qs[2] + 1.5 * iqr)
  out
}

sig_symbol <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

make_long_data <- function(status_col, status_levels) {
  state_df %>%
    select(Run, all_of(colnames(state_fraction)), status = all_of(status_col)) %>%
    filter(status %in% status_levels) %>%
    mutate(status = factor(status, levels = status_levels)) %>%
    pivot_longer(cols = all_of(colnames(state_fraction)), names_to = "CellType", values_to = "Abundance") %>%
    mutate(Abundance = as.numeric(Abundance))
}

filter_iqr_by_cell_status <- function(dat) {
  dat %>%
    group_by(CellType, status) %>%
    mutate(IsOutlier = mark_iqr_outlier(Abundance)) %>%
    ungroup() %>%
    filter(!IsOutlier)
}

calc_pair_stats <- function(dat, comparisons) {
  pieces <- list()
  k <- 1L
  for (cell_i in sort(unique(dat$CellType))) {
    sub <- dat[dat$CellType == cell_i, , drop = FALSE]
    for (cmp in comparisons) {
      v1 <- sub$Abundance[sub$status == cmp[1]]
      v2 <- sub$Abundance[sub$status == cmp[2]]
      v1 <- v1[is.finite(v1)]
      v2 <- v2[is.finite(v2)]
      p <- if (length(v1) > 0 && length(v2) > 0) {
        suppressWarnings(wilcox.test(v1, v2)$p.value)
      } else {
        NA_real_
      }
      pieces[[k]] <- data.frame(
        CellType = cell_i,
        group1 = cmp[1],
        group2 = cmp[2],
        n1 = length(v1),
        n2 = length(v2),
        mean1 = mean(v1, na.rm = TRUE),
        mean2 = mean(v2, na.rm = TRUE),
        median1 = median(v1, na.rm = TRUE),
        median2 = median(v2, na.rm = TRUE),
        log2FC_group2_vs_group1 = log2((mean(v2, na.rm = TRUE) + 1e-12) / (mean(v1, na.rm = TRUE) + 1e-12)),
        p = p,
        p.signif = sig_symbol(p),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  bind_rows(pieces) %>%
    group_by(group1, group2) %>%
    mutate(FDR = p.adjust(p, method = "BH")) %>%
    ungroup()
}

make_half_violin_data <- function(dat, x_map, width = 0.25) {
  dat <- dat %>% filter(!is.na(Abundance), !is.na(status), !is.na(CellType))
  pieces <- list()
  k <- 1L
  for (status_i in names(x_map)) {
    vals <- dat$Abundance[as.character(dat$status) == status_i]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L || length(unique(vals)) < 2L) next
    dens <- density(vals, n = 128, from = min(vals), to = max(vals), na.rm = TRUE)
    dens_scaled <- dens$y / max(dens$y, na.rm = TRUE) * width
    x_base <- unname(x_map[status_i]) + 0.20
    pieces[[k]] <- data.frame(
      status = factor(status_i, levels = names(x_map)),
      x = c(x_base + dens_scaled, rep(x_base, length(dens$x))),
      y = c(dens$x, rev(dens$x)),
      group_id = status_i,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  bind_rows(pieces)
}

make_cell_panel <- function(dat, stat_df, cell, comparisons, x_map, show_y = FALSE) {
  tmp <- dat %>%
    filter(CellType == cell) %>%
    mutate(
      x = unname(x_map[as.character(status)]),
      x_box = x - 0.16
    )
  half_df <- make_half_violin_data(tmp, x_map = x_map)
  y_max <- max(tmp$Abundance, na.rm = TRUE)
  y_min <- min(tmp$Abundance, na.rm = TRUE)
  y_range <- ifelse(y_max > y_min, y_max - y_min, max(y_max, 1e-6))
  y_pad <- y_range * 0.22
  stat_plot <- stat_df %>%
    filter(CellType == cell) %>%
    mutate(
      x = unname(x_map[group1]),
      xend = unname(x_map[group2]),
      y = y_max + seq_along(p) * y_range * 0.10,
      yend = y,
      label_y = y + y_range * 0.035
    )

  ggplot(tmp, aes(x = x_box, y = Abundance)) +
    geom_polygon(
      data = half_df,
      aes(x = x, y = y, group = group_id, fill = status),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.75,
      alpha = 0.9
    ) +
    stat_boxplot(
      aes(group = status, colour = status),
      geom = "errorbar",
      width = 0.16,
      linewidth = 0.7
    ) +
    geom_boxplot(
      aes(group = status, colour = status, fill = status),
      width = 0.32,
      outlier.shape = NA,
      linewidth = 0.7,
      alpha = 0.95
    ) +
    stat_summary(
      aes(group = status),
      fun = "median",
      geom = "point",
      shape = 16,
      size = 1.8,
      color = "white"
    ) +
    geom_segment(
      data = stat_plot,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.65
    ) +
    geom_text(
      data = stat_plot,
      aes(x = (x + xend) / 2, y = label_y, label = p.signif),
      inherit.aes = FALSE,
      size = 5.5,
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = unname(x_map),
      labels = names(x_map),
      limits = c(min(unname(x_map)) - 0.55, max(unname(x_map)) + 0.65),
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = c(
      "Died" = "#F28E2BFF",
      "Survived" = "#76B7B2FF",
      "Died-IRIS" = "#F28E2BFF",
      "Survived-IRIS" = "#76B7B2FF",
      "Early ART" = "#00A087FF",
      "Deferred ART" = "#EFC000FF"
    )) +
    scale_color_manual(values = setNames(rep("black", length(x_map)), names(x_map))) +
    labs(x = NULL, y = if (show_y) "Cell Ratio" else NULL, title = cell) +
    coord_cartesian(ylim = c(y_min, y_max + y_pad + nrow(stat_plot) * y_range * 0.08), clip = "off") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", color = "#8B1A1A", hjust = 0),
      axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.y = element_text(size = 12, color = "#8B1A1A", face = "bold"),
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.75, linetype = "dashed"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.95),
      aspect.ratio = 1,
      plot.margin = margin(4, 8, 12, 8)
    )
}

make_multi_panel_plot <- function(dat, stat_df, selected_cells, comparisons, status_levels, output_prefix,
                                  nrow = 1, width = 12, height = 3.6) {
  selected_cells <- selected_cells[selected_cells %in% unique(dat$CellType)]
  dat <- dat %>% filter(CellType %in% selected_cells, status %in% status_levels) %>%
    mutate(status = factor(as.character(status), levels = status_levels))
  stat_df <- stat_df %>% filter(CellType %in% selected_cells)
  x_map <- setNames(seq_along(status_levels), status_levels)
  plot_list <- Map(
    function(cell, idx) make_cell_panel(dat, stat_df, cell, comparisons, x_map, show_y = idx == 1L),
    selected_cells,
    seq_along(selected_cells)
  )
  p <- cowplot::plot_grid(
    plotlist = plot_list,
    nrow = nrow,
    align = "hv",
    labels = c("A", rep("", length(plot_list) - 1L)),
    label_size = 18,
    label_fontface = "bold",
    hjust = -0.2,
    vjust = 1.4
  )
  ggsave(file.path(out_dir, paste0(output_prefix, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(out_dir, paste0(output_prefix, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
  p
}

message("Preparing outcome comparisons...")
outcome_levels <- c("Died", "Survived", "Died-IRIS", "Survived-IRIS")
outcome_comparisons <- list(c("Died", "Survived"), c("Died-IRIS", "Survived-IRIS"))
outcome_long <- make_long_data("outcome", outcome_levels)
outcome_filtered <- filter_iqr_by_cell_status(outcome_long)
outcome_stats <- calc_pair_stats(outcome_filtered, outcome_comparisons)

write.table(outcome_long, file.path(out_dir, "PRJNA683803_outcome_all_cells_long.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(outcome_filtered, file.path(out_dir, "PRJNA683803_outcome_all_cells_long_no_outliers.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(outcome_stats, file.path(out_dir, "PRJNA683803_outcome_all_celltype_wilcox_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

died_survived_selected <- c("MDSC", "Neutrophil", "Th1", "Tc")
died_survived_stats <- outcome_stats %>%
  filter(group1 == "Died", group2 == "Survived")
write.table(died_survived_stats %>% filter(CellType %in% died_survived_selected),
            file.path(out_dir, "PRJNA683803_Died_vs_Survived_selected_cell_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
make_multi_panel_plot(
  outcome_filtered,
  died_survived_stats,
  selected_cells = died_survived_selected,
  comparisons = list(c("Died", "Survived")),
  status_levels = c("Died", "Survived"),
  output_prefix = "Figure6_PRJNA683803_Died_vs_Survived_selected_cells",
  nrow = 1,
  width = 12.5,
  height = 3.2
)

iris_selected <- c("MDSC", "Neutrophil", "Th1", "Tc")
iris_stats <- outcome_stats %>%
  filter(group1 == "Died-IRIS", group2 == "Survived-IRIS")
write.table(iris_stats %>% filter(CellType %in% iris_selected),
            file.path(out_dir, "PRJNA683803_DiedIRIS_vs_SurvivedIRIS_selected_cell_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
make_multi_panel_plot(
  outcome_filtered,
  iris_stats,
  selected_cells = iris_selected,
  comparisons = list(c("Died-IRIS", "Survived-IRIS")),
  status_levels = c("Died-IRIS", "Survived-IRIS"),
  output_prefix = "Figure6_PRJNA683803_DiedIRIS_vs_SurvivedIRIS_selected_cells",
  nrow = 1,
  width = 12.5,
  height = 3.2
)

message("Preparing ART timing comparison...")
treatment_levels <- c("Early ART", "Deferred ART")
treatment_comparisons <- list(c("Early ART", "Deferred ART"))
treatment_long <- make_long_data("treatment", treatment_levels)
treatment_filtered <- filter_iqr_by_cell_status(treatment_long)
treatment_stats <- calc_pair_stats(treatment_filtered, treatment_comparisons)

write.table(treatment_long, file.path(out_dir, "PRJNA683803_treatment_all_cells_long.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(treatment_filtered, file.path(out_dir, "PRJNA683803_treatment_all_cells_long_no_outliers.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(treatment_stats, file.path(out_dir, "PRJNA683803_treatment_all_celltype_wilcox_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

treatment_selected <- c("MDSC", "CD4Tn", "Th1", "gdT")
write.table(treatment_stats %>% filter(CellType %in% treatment_selected),
            file.path(out_dir, "PRJNA683803_EarlyART_vs_DeferredART_selected_cell_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
make_multi_panel_plot(
  treatment_filtered,
  treatment_stats,
  selected_cells = treatment_selected,
  comparisons = treatment_comparisons,
  status_levels = treatment_levels,
  output_prefix = "Figure6_PRJNA683803_ART_Early_vs_Deferred_selected_cells",
  nrow = 1,
  width = 12.5,
  height = 3.2
)

treatment_supplement_selected <- c("CD8Tn", "Th1/17", "Tr1", "cNK")
write.table(treatment_stats %>% filter(CellType %in% treatment_supplement_selected),
            file.path(out_dir, "PRJNA683803_EarlyART_vs_DeferredART_supplement_selected_cell_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
make_multi_panel_plot(
  treatment_filtered,
  treatment_stats,
  selected_cells = treatment_supplement_selected,
  comparisons = treatment_comparisons,
  status_levels = treatment_levels,
  output_prefix = "Figure6_PRJNA683803_ART_Early_vs_Deferred_supplement_selected_cells",
  nrow = 1,
  width = 12.5,
  height = 3.2
)

run_info <- data.frame(
  Item = c("Expression", "Metadata", "Reference", "MarkerFile", "Samples",
           "OutcomeCounts", "TreatmentCounts", "GenesUsed", "MarkerTokens",
           "MarkerOverlap", "HierarchyMode", "InferenceMethod", "AddUnknown",
           "Cores", "PackageVersion"),
  Value = c(
    exp_file,
    metadata_file,
    reference_file,
    immune_gene_file,
    nrow(metadata),
    paste(names(table(metadata$outcome)), as.integer(table(metadata$outcome)), sep = "=", collapse = "; "),
    paste(names(table(metadata$treatment)), as.integer(table(metadata$treatment)), sep = "=", collapse = "; "),
    genes_used,
    marker_tokens_n,
    marker_overlap_n,
    "tcell",
    "vb",
    "FALSE",
    n_cores,
    as.character(packageVersion("ImmuCellAI2.0"))
  )
)
write.table(run_info, file.path(out_dir, "ImmuCellAI2_PRJNA683803_run_info.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done. PRJNA683803 AIDS results written to:\n", out_dir, "\n", sep = "")
print(run_info)
