# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

root <- file.path(fig3_work_dir(), "Fig2/Newman_Monaco_7tools_validation")
out_dir <- file.path(root, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("ggplot2", "RColorBrewer")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing package(s): ", paste(missing, collapse = ", "))

library(ggplot2)

inputs <- data.frame(
  Dataset = c("Newman PBMC", "Monaco PBMC"),
  File = c(
    file.path(root, "RNA_Sieve_newman_pbmcs", "seven_tools_compare",
              "all_methods_per_celltype_metrics.txt"),
    file.path(root, "RNA_Sieve_monaco_pbmcs", "seven_tools_compare",
              "all_methods_per_celltype_metrics.txt")
  ),
  stringsAsFactors = FALSE
)

method_labels <- c(
  ImmuCellAI2_flat_VB = "ImmuCellAI 2.0",
  ImmuCellAI2_tcell_VB_UNKNOWN = "ImmuCellAI 2.0",
  BayesPrism_first_state_chain600_burn500 = "BayesPrism",
  CIBERSORT_default_LM22 = "CIBERSORT",
  CITMIC_native = "CITMIC",
  DWLS_native_ref_top1000 = "DWLS",
  ImmuCellAI_native = "ImmuCellAI",
  MuSiC_native_music.basic = "MuSiC"
)

method_order <- c(
  "ImmuCellAI 2.0", "BayesPrism", "CIBERSORT", "CITMIC",
  "DWLS", "ImmuCellAI", "MuSiC"
)

cell_labels <- c(
  B = "B cell",
  CD4T = "CD4 T cell",
  CD8T = "CD8 T cell",
  Monocytes = "Monocyte",
  Neutrophils = "Neutrophil",
  NK = "NK cell"
)

read_one <- function(dataset, file) {
  x <- read.delim(file, check.names = FALSE)
  x$Dataset <- dataset
  x$method <- unname(method_labels[x$Method])
  x$cellType <- ifelse(
    x$CellType %in% names(cell_labels),
    unname(cell_labels[x$CellType]),
    x$CellType
  )
  x$correlation <- suppressWarnings(as.numeric(x$Pearson))
  x$status <- "valid"
  x$status[!is.finite(x$correlation) & trimws(x$PredictedStates) == ""] <- "not available"
  x$status[!is.finite(x$correlation) & trimws(x$PredictedStates) != ""] <- "constant prediction"
  x$correlation[!is.finite(x$correlation)] <- NA_real_
  x
}

plot_data <- do.call(rbind, Map(read_one, inputs$Dataset, inputs$File))
write.csv(plot_data, file.path(out_dir, "Newman_Monaco_7tools_pie_plot_data.csv"),
          row.names = FALSE, quote = FALSE)

theme_blue <- theme(
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"),
  axis.text.y = element_text(size = 12, colour = "black"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.65),
  strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 0.8),
  strip.text = element_text(size = 18, face = "bold", colour = "black"),
  legend.title = element_text(size = 12, colour = "black"),
  legend.text = element_text(size = 11, colour = "black"),
  legend.key = element_blank(),
  plot.margin = margin(8, 12, 8, 8)
)

make_circle <- function(cx, cy, r = 0.32, n = 100L) {
  a <- seq(0, 2 * pi, length.out = n + 1L)
  data.frame(x = cx + r * cos(a), y = cy + r * sin(a))
}

make_wedge <- function(cx, cy, correlation, r = 0.32, n = 100L) {
  frac <- min(abs(correlation), 1)
  a <- seq(pi / 2, pi / 2 - 2 * pi * frac,
           length.out = max(3L, ceiling(n * frac) + 1L))
  data.frame(
    x = c(cx, cx + r * cos(a), cx),
    y = c(cy, cy + r * sin(a), cy)
  )
}

draw_dataset <- function(dataset) {
  d <- plot_data[plot_data$Dataset == dataset & !is.na(plot_data$method), , drop = FALSE]
  cell_order <- sort(unique(d$cellType))
  d$cellType <- factor(d$cellType, levels = cell_order)
  d$method <- factor(d$method, levels = rev(method_order))
  d$x_id <- as.numeric(d$cellType)
  d$y_id <- as.numeric(d$method)

  valid <- d[is.finite(d$correlation), , drop = FALSE]
  constant <- d[d$status == "constant prediction", , drop = FALSE]
  unavailable <- d[d$status == "not available", , drop = FALSE]

  circles <- do.call(rbind, lapply(seq_len(nrow(valid)), function(i) {
    z <- make_circle(valid$x_id[i], valid$y_id[i])
    z$id <- i
    z$Dataset <- dataset
    z
  }))
  wedges <- do.call(rbind, lapply(seq_len(nrow(valid)), function(i) {
    z <- make_wedge(valid$x_id[i], valid$y_id[i], valid$correlation[i])
    z$id <- i
    z$correlation <- valid$correlation[i]
    z$Dataset <- dataset
    z
  }))

  p <- ggplot() +
    geom_polygon(
      data = circles,
      aes(x = x, y = y, group = id),
      fill = "white",
      color = "#1F1F1F",
      linewidth = 0.22
    ) +
    geom_polygon(
      data = wedges,
      aes(x = x, y = y, group = id, fill = correlation),
      color = "#1F1F1F",
      linewidth = 0.18
    ) +
    geom_point(
      data = constant,
      aes(x = x_id, y = y_id),
      inherit.aes = FALSE,
      shape = 4,
      size = 2.2,
      color = "black",
      stroke = 0.55
    ) +
    geom_text(
      data = unavailable,
      aes(x = x_id, y = y_id),
      inherit.aes = FALSE,
      label = "/",
      size = 4,
      fontface = "bold",
      color = "black"
    ) +
    facet_grid(. ~ Dataset) +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      limits = c(-1, 1),
      name = "Pearson"
    ) +
    scale_x_continuous(
      breaks = seq_along(cell_order),
      labels = cell_order,
      limits = c(0.35, length(cell_order) + 0.65),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(rev(method_order)),
      labels = rev(method_order),
      limits = c(0.35, length(method_order) + 0.65),
      expand = c(0, 0)
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    theme_blue

  safe_name <- gsub(" ", "_", dataset)
  width <- max(7.2, 1.0 * length(cell_order) + 2.8)
  ggsave(file.path(out_dir, paste0(safe_name, "_7tools_per_celltype_pie.pdf")),
         p, width = width, height = 6.2)
  ggsave(file.path(out_dir, paste0(safe_name, "_7tools_per_celltype_pie.png")),
         p, width = width, height = 6.2, dpi = 300)
  p
}

plots <- lapply(inputs$Dataset, draw_dataset)
names(plots) <- inputs$Dataset

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package patchwork is required for the combined plot.")
}

compact_theme <- theme(
  strip.text = element_text(size = 8.5, face = "bold", colour = "black"),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.5, colour = "black"),
  axis.text.y = element_text(size = 5.2, colour = "black"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
  strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 0.25),
  axis.ticks = element_line(linewidth = 0.15),
  legend.title = element_text(size = 6.5, colour = "black"),
  legend.text = element_text(size = 6, colour = "black"),
  legend.key.height = grid::unit(2.5, "mm"),
  legend.key.width = grid::unit(4.5, "mm"),
  plot.margin = margin(2, 2, 2, 2, "mm")
)

combined_plot <- patchwork::wrap_plots(
  plots[["Newman PBMC"]] + compact_theme,
  plots[["Monaco PBMC"]] + compact_theme,
  ncol = 2,
  guides = "collect",
  widths = c(1.08, 1)
) & theme(legend.position = "bottom")

ggsave(
  file.path(out_dir, "Newman_Monaco_7tools_per_celltype_pies_combined.pdf"),
  combined_plot,
  width = 4.96,
  height = 3.1,
  units = "in",
  bg = "white",
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "Newman_Monaco_7tools_per_celltype_pies_combined.png"),
  combined_plot,
  width = 4.96,
  height = 3.1,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

message("Plots written to: ", out_dir)
