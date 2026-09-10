# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

work_dir <- fig3_work_dir()
fig2_dir <- file.path(work_dir, "Fig2")
out_dir <- file.path(fig2_dir, "GSE107011_Newman_Monaco_runtime_combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gse_file <- file.path(
  work_dir, "combined_4datasets_7tools_pie_plots_final_ordered",
  "GSE107011_7tools_Figure2A_pie_plot_data.csv"
)
pbmc_file <- file.path(
  fig2_dir, "Newman_Monaco_7tools_validation", "plots",
  "Newman_Monaco_7tools_pie_plot_data.csv"
)
runtime_file <- file.path(
  fig2_dir, "runtime_compare",
  "runtime_benchmark_1sample_vs_50samples_summary.txt"
)

missing_files <- c(gse_file, pbmc_file, runtime_file)
missing_files <- missing_files[!file.exists(missing_files)]
if (length(missing_files)) stop("Missing input files:\n", paste(missing_files, collapse = "\n"))

required <- c("ggplot2", "RColorBrewer", "patchwork")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))

library(ggplot2)
library(patchwork)

method_order <- c(
  "ImmuCellAI 2.0", "BayesPrism", "CIBERSORT", "CITMIC",
  "DWLS", "ImmuCellAI", "MuSiC"
)

method_map <- c(
  "ImmuCellAI 2.0" = "ImmuCellAI 2.0",
  "ImmuCellAI 2.0" = "ImmuCellAI 2.0",
  BayesPrism = "BayesPrism",
  CIBERSORT = "CIBERSORT",
  CITMIC = "CITMIC",
  DWLS = "DWLS",
  ImmuCellAI = "ImmuCellAI",
  MuSiC = "MuSiC"
)

gse <- read.csv(gse_file, check.names = FALSE)
gse$Dataset <- "GSE107011"
gse$method <- unname(method_map[as.character(gse$method)])
gse$cellType <- as.character(gse$cellTypeLabel)
gse$correlation <- suppressWarnings(as.numeric(gse$correlation))
gse$status <- as.character(gse$status)

pbmc <- read.csv(pbmc_file, check.names = FALSE)
pbmc$method <- unname(method_map[as.character(pbmc$method)])
pbmc$cellType <- as.character(pbmc$cellType)
pbmc$correlation <- suppressWarnings(as.numeric(pbmc$correlation))
pbmc$status <- gsub(" ", "_", as.character(pbmc$status), fixed = TRUE)

pie_data <- rbind(
  gse[, c("Dataset", "method", "cellType", "correlation", "status")],
  pbmc[, c("Dataset", "method", "cellType", "correlation", "status")]
)
pie_data <- pie_data[pie_data$method %in% method_order, , drop = FALSE]

write.csv(
  pie_data,
  file.path(out_dir, "combined_pie_plot_data.csv"),
  row.names = FALSE,
  quote = FALSE
)

make_circle <- function(cx, cy, r = 0.31, n = 80L) {
  angle <- seq(0, 2 * pi, length.out = n + 1L)
  data.frame(x = cx + r * cos(angle), y = cy + r * sin(angle))
}

make_wedge <- function(cx, cy, value, r = 0.31, n = 80L) {
  fraction <- min(max(abs(value), 0), 1)
  if (!is.finite(fraction) || fraction <= 0) return(NULL)
  angle <- seq(
    pi / 2,
    pi / 2 - 2 * pi * fraction,
    length.out = max(3L, ceiling(n * fraction) + 1L)
  )
  data.frame(
    x = c(cx, cx + r * cos(angle), cx),
    y = c(cy, cy + r * sin(angle), cy)
  )
}

make_pie_polygons <- function(valid, radius = 0.31) {
  circles <- vector("list", nrow(valid))
  wedges <- vector("list", nrow(valid))
  for (i in seq_len(nrow(valid))) {
    circle <- make_circle(valid$x_id[i], valid$y_id[i], radius)
    circle$id <- i
    circles[[i]] <- circle
    wedge <- make_wedge(valid$x_id[i], valid$y_id[i], valid$correlation[i], radius)
    if (!is.null(wedge)) {
      wedge$id <- i
      wedge$correlation <- valid$correlation[i]
    }
    wedges[[i]] <- wedge
  }
  list(
    circles = if (length(circles)) do.call(rbind, circles) else data.frame(),
    wedges = if (length(Filter(Negate(is.null), wedges))) {
      do.call(rbind, Filter(Negate(is.null), wedges))
    } else {
      data.frame()
    }
  )
}

draw_pie_panel <- function(dataset, top_panel = FALSE) {
  d <- pie_data[pie_data$Dataset == dataset, , drop = FALSE]
  d$method <- factor(d$method, levels = rev(method_order))
  cell_order <- sort(unique(as.character(d$cellType)))
  d$cellType <- factor(d$cellType, levels = cell_order)
  d$x_id <- as.numeric(d$cellType)
  d$y_id <- as.numeric(d$method)

  valid <- d[is.finite(d$correlation), , drop = FALSE]
  constant <- d[d$status == "constant_prediction", , drop = FALSE]
  unavailable <- d[d$status == "not_available", , drop = FALSE]
  polygons <- make_pie_polygons(valid, radius = 0.31)
  polygons$circles$Dataset <- dataset
  polygons$wedges$Dataset <- dataset

  x_size <- if (top_panel) 4.7 else 5.0
  y_size <- if (top_panel) 5.7 else 5.5
  title_size <- if (top_panel) 9.2 else 8.2

  ggplot() +
    geom_polygon(
      data = polygons$circles,
      aes(x = x, y = y, group = id),
      inherit.aes = FALSE,
      fill = "white",
      color = "#202020",
      linewidth = 0.16
    ) +
    geom_polygon(
      data = polygons$wedges,
      aes(x = x, y = y, group = id, fill = correlation),
      inherit.aes = FALSE,
      color = "#202020",
      linewidth = 0.13
    ) +
    geom_point(
      data = constant,
      aes(x = x_id, y = y_id),
      inherit.aes = FALSE,
      shape = 4,
      size = 0.85,
      stroke = 0.45,
      color = "black"
    ) +
    geom_text(
      data = unavailable,
      aes(x = x_id, y = y_id),
      inherit.aes = FALSE,
      label = "/",
      size = 2.2,
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
      expand = expansion(add = 0.48)
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(d$method)),
      labels = levels(d$method),
      expand = expansion(add = 0.45)
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 7) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.20),
      strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 0.20),
      strip.text = element_text(size = title_size, face = "bold", colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = x_size, colour = "black"),
      axis.text.y = element_text(size = y_size, colour = "black"),
      axis.ticks = element_line(linewidth = 0.12),
      legend.title = element_text(size = 6.2, colour = "black"),
      legend.text = element_text(size = 5.8, colour = "black"),
      legend.key.height = grid::unit(2.4, "mm"),
      legend.key.width = grid::unit(4.5, "mm"),
      plot.margin = margin(1.4, 1.8, 1.4, 1.8, "mm")
    )
}

runtime <- read.delim(runtime_file, check.names = FALSE)
runtime_map <- c(
  ImmuCellAI2_tcell_VB = "ImmuCellAI 2.0",
  BayesPrism_first_state_chain600_burn500 = "BayesPrism",
  CIBERSORT_nuSVR_ref_top1000 = "CIBERSORT",
  CITMIC_native = "CITMIC",
  DWLS_weighted_lm = "DWLS",
  ImmuCellAI_native = "ImmuCellAI",
  MuSiC_basic = "MuSiC"
)
runtime$MethodShort <- unname(runtime_map[as.character(runtime$Method)])
runtime <- runtime[runtime$Status == "OK" & runtime$MethodShort %in% method_order, , drop = FALSE]
runtime$MethodShort <- factor(runtime$MethodShort, levels = rev(method_order))
runtime$ScenarioLabel <- factor(
  runtime$Scenario,
  levels = c("single_sample_1core", "fifty_samples_8cores"),
  labels = c("1 sample / 1 core", "50 samples / 8 cores")
)
runtime$LogTime <- log10(suppressWarnings(as.numeric(runtime$SecondsPerSample)) + 1)
runtime$Panel <- "Runtime efficiency"

write.table(
  runtime,
  file.path(out_dir, "runtime_plot_data.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_gse <- draw_pie_panel("GSE107011", top_panel = TRUE) +
  theme(legend.position = "right")

p_newman <- draw_pie_panel("Newman PBMC") +
  theme(legend.position = "none")

p_monaco <- draw_pie_panel("Monaco PBMC") +
  theme(legend.position = "none")

p_runtime <- ggplot(runtime, aes(x = MethodShort, y = LogTime, fill = ScenarioLabel)) +
  geom_col(
    position = position_dodge(width = 0.68),
    width = 0.60,
    color = "#202020",
    linewidth = 0.16
  ) +
  coord_flip() +
  facet_grid(. ~ Panel) +
  scale_fill_manual(values = c("1 sample / 1 core" = "#379DA5", "50 samples / 8 cores" = "#F66463")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "log10(seconds per sample + 1)", fill = NULL) +
  theme_bw(base_size = 7) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.25, linetype = "dashed"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.20),
    strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 0.20),
    strip.text = element_text(size = 8.2, face = "bold", colour = "black"),
    axis.text = element_text(size = 5.5, colour = "black"),
    axis.title.x = element_text(size = 5.8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.12),
    legend.position = "bottom",
    legend.text = element_text(size = 5.3, colour = "black"),
    legend.key.height = grid::unit(2.3, "mm"),
    legend.key.width = grid::unit(4.2, "mm"),
    plot.margin = margin(1.4, 1.8, 1.4, 1.8, "mm")
  )

bottom_row <- p_newman + p_monaco + p_runtime +
  plot_layout(ncol = 3, widths = c(1.05, 1, 0.92))

combined_plot <- p_gse / bottom_row +
  plot_layout(heights = c(1.02, 1))

out_pdf <- file.path(out_dir, "Figure2_GSE107011_Newman_Monaco_runtime_combined.pdf")
out_png <- file.path(out_dir, "Figure2_GSE107011_Newman_Monaco_runtime_combined.png")

ggsave(out_pdf, combined_plot, width = 11.69, height = 8.27, units = "in", bg = "white", limitsize = FALSE)
ggsave(out_png, combined_plot, width = 11.69, height = 8.27, units = "in", dpi = 500, bg = "white", limitsize = FALSE)

message("Saved combined PDF: ", out_pdf)
message("Saved combined PNG: ", out_png)
