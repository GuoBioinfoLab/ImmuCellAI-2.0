options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
})

fig6_dir <- "<LOCAL_R_ROOT>/Fig6"
out_dir <- file.path(fig6_dir, "ImmuCellAI2_active_latent_TB")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

exp_file <- file.path(fig6_dir, "exp.txt")
metadata_file <- file.path(fig6_dir, "PRJNA_combined.csv")
reference_file <- "<LOCAL_R_ROOT>/reference_53celltypesTPM20260518.txt"
immune_gene_file <- "<LOCAL_R_ROOT>/MarkerUsedDeconvolution.txt"
package_file <- "<LOCAL_R_ROOT>/ImmuCellAI2.0_0.1.7.tar.gz"

stable_tmp <- file.path(out_dir, "Rtmp_immucellai2_tb")
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

message("Reading TB sample metadata...")
metadata <- fread(metadata_file, data.table = FALSE, check.names = FALSE) %>%
  filter(tb_status %in% c("Active TB", "Latent TB")) %>%
  distinct(Run, .keep_all = TRUE)
metadata$tb_status <- factor(metadata$tb_status, levels = c("Active TB", "Latent TB"))

message("Checking expression columns...")
exp_header <- names(fread(exp_file, nrows = 0, data.table = FALSE, check.names = FALSE))
matched_samples <- metadata$Run[metadata$Run %in% exp_header]
if (length(matched_samples) == 0) {
  stop("No metadata samples are present in exp.txt.")
}
metadata <- metadata[match(matched_samples, metadata$Run), , drop = FALSE]

message("Reading selected expression columns: ", length(matched_samples), " samples...")
exp_df <- fread(exp_file, data.table = FALSE, check.names = FALSE,
                select = c("V1", matched_samples))
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

state_fraction_file <- file.path(out_dir, "ImmuCellAI2_TB_state_fraction_sample_by_celltype.txt")
if (file.exists(state_fraction_file)) {
  message("Using existing ImmuCellAI2 result: ", state_fraction_file)
  state_fraction <- read.table(state_fraction_file, header = TRUE, sep = "\t", row.names = 1,
                               check.names = FALSE, comment.char = "")
  state_fraction <- as.matrix(state_fraction)
  state_fraction <- state_fraction[matched_samples, , drop = FALSE]
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
  state_fraction <- state_fraction[matched_samples, , drop = FALSE]
  write.table(state_fraction, state_fraction_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(t(state_fraction), file.path(out_dir, "ImmuCellAI2_TB_state_fraction_celltype_by_sample.txt"),
              sep = "\t", quote = FALSE, col.names = NA)
}

selected_map <- c(cDC1 = "cDC1", cMo = "cMo", MDSC = "MDSC", Tc = "Tc")
missing_cells <- setdiff(unname(selected_map), colnames(state_fraction))
if (length(missing_cells) > 0) {
  stop("Missing selected ImmuCellAI2 cells: ", paste(missing_cells, collapse = ", "))
}

plot_df <- as.data.frame(state_fraction[, unname(selected_map), drop = FALSE]) %>%
  tibble::rownames_to_column("Run")
colnames(plot_df)[match(unname(selected_map), colnames(plot_df))] <- names(selected_map)
plot_df <- plot_df %>%
  left_join(metadata %>% select(Run, BioProject, Disease, tb_status), by = "Run") %>%
  pivot_longer(cols = all_of(names(selected_map)), names_to = "CellType", values_to = "Abundance") %>%
  mutate(
    status = factor(tb_status, levels = c("Active TB", "Latent TB")),
    CellType = factor(CellType, levels = c("cDC1", "Tc", "cMo", "MDSC")),
    Abundance = as.numeric(Abundance)
  )

write.table(metadata, file.path(out_dir, "ImmuCellAI2_TB_matched_metadata.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(plot_df, file.path(out_dir, "ImmuCellAI2_TB_Active_vs_Latent_selected_cells_long.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

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

remove_plot_outlier_cells <- character(0)
plot_df_for_plot <- plot_df %>%
  group_by(CellType, status) %>%
  mutate(
    IsPlotOutlier = if (as.character(first(CellType)) %in% remove_plot_outlier_cells) {
      mark_iqr_outlier(Abundance)
    } else {
      rep(FALSE, dplyr::n())
    }
  ) %>%
  ungroup()

plot_df_no_outliers <- plot_df_for_plot %>%
  filter(!IsPlotOutlier)

write.table(plot_df_for_plot, file.path(out_dir, "ImmuCellAI2_TB_Active_vs_Latent_selected_cells_long_with_outlier_flag.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(plot_df_no_outliers, file.path(out_dir, "ImmuCellAI2_TB_Active_vs_Latent_selected_cells_long_for_plot.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

stats <- plot_df_no_outliers %>%
  group_by(CellType) %>%
  summarise(
    N_Active = sum(status == "Active TB" & !is.na(Abundance)),
    N_Latent = sum(status == "Latent TB" & !is.na(Abundance)),
    Mean_Active = mean(Abundance[status == "Active TB"], na.rm = TRUE),
    Mean_Latent = mean(Abundance[status == "Latent TB"], na.rm = TRUE),
    P_Wilcox = suppressWarnings(wilcox.test(Abundance ~ status)$p.value),
    .groups = "drop"
  )
write.table(stats, file.path(out_dir, "ImmuCellAI2_TB_Active_vs_Latent_selected_cells_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sig_symbol <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

make_half_violin_data <- function(dat, width = 0.25) {
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
    x_base <- ifelse(status_i == "Active TB", 1.20, 2.20)
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

make_tb_panel <- function(cell, show_y = FALSE) {
  tmp <- plot_df_no_outliers %>%
    filter(CellType == cell) %>%
    mutate(
      x = ifelse(status == "Active TB", 1, 2),
      x_box = x - 0.16
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
      limits = c(0.48, 2.62),
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

panel_cells <- c("cDC1", "Tc", "cMo", "MDSC")
plot_list <- Map(make_tb_panel, panel_cells, c(TRUE, FALSE, FALSE, FALSE))
p_combined <- cowplot::plot_grid(plotlist = plot_list, nrow = 1, align = "hv",
                                 labels = c("A", "", "", ""), label_size = 18,
                                 label_fontface = "bold", hjust = -0.2, vjust = 1.4)

pdf_file <- file.path(out_dir, "Figure6A_ImmuCellAI2_ActiveTB_vs_LatentTB_selected_cells.pdf")
png_file <- file.path(out_dir, "Figure6A_ImmuCellAI2_ActiveTB_vs_LatentTB_selected_cells.png")
ggsave(pdf_file, p_combined, width = 12.5, height = 3.2, bg = "white")
ggsave(png_file, p_combined, width = 12.5, height = 3.2, dpi = 600, bg = "white")

run_info <- data.frame(
  Item = c("Expression", "Metadata", "Reference", "MarkerFile", "Samples",
           "ActiveTB", "LatentTB", "GenesUsed", "MarkerTokens", "MarkerOverlap",
           "HierarchyMode", "InferenceMethod", "AddUnknown", "Cores", "PackageVersion"),
  Value = c(exp_file, metadata_file, reference_file, immune_gene_file, nrow(metadata),
            sum(metadata$tb_status == "Active TB"), sum(metadata$tb_status == "Latent TB"),
            nrow(bulk), length(marker_tokens), length(marker_genes),
            "tcell", "vb", "FALSE", n_cores, as.character(packageVersion("ImmuCellAI2.0")))
)
write.table(run_info, file.path(out_dir, "ImmuCellAI2_TB_run_info.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done. Results written to:\n", out_dir, "\n", sep = "")
print(stats)
