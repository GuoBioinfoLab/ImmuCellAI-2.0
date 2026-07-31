suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scRNAtoolVis)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    paste(
      "Usage: Rscript benchmarks/plot_per_celltype_correlations.R",
      "combined_metrics.tsv output.pdf [pie_width]"
    )
  )
}

input_file <- args[1]
output_file <- args[2]
pie_width <- if (length(args) >= 3L) as.numeric(args[3]) else 0.65
if (!is.finite(pie_width) || pie_width <= 0) pie_width <- 0.65

dat <- fread(input_file, data.table = FALSE, check.names = FALSE)
required <- c("Dataset", "Method", "TargetCell", "Pearson", "Status")
missing_columns <- setdiff(required, names(dat))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

names(dat)[names(dat) == "TargetCell"] <- "cellType"
names(dat)[names(dat) == "Method"] <- "method"
names(dat)[names(dat) == "Pearson"] <- "correlation"

method_levels <- c(
  "ImmuCellAI 2.0",
  sort(setdiff(unique(dat$method), "ImmuCellAI 2.0"))
)
dat$method <- factor(dat$method, levels = rev(method_levels))
dat$cellType <- factor(dat$cellType, levels = sort(unique(dat$cellType)))
dat$Dataset <- factor(dat$Dataset, levels = unique(dat$Dataset))

valid <- dat[dat$Status == "valid" & is.finite(dat$correlation), , drop = FALSE]
constant <- dat[dat$Status == "constant_prediction", , drop = FALSE]
unavailable <- dat[dat$Status == "unavailable", , drop = FALSE]
insufficient <- dat[dat$Status == "insufficient_samples", , drop = FALSE]

spectral <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
p <- ggplot(dat, aes(x = cellType, y = method, fill = correlation)) +
  scRNAtoolVis::geom_jjpie(
    data = valid,
    aes(piefill = correlation),
    width = pie_width,
    color = "black",
    linewidth = 0.45
  ) +
  geom_point(
    data = constant,
    shape = 4,
    size = 2,
    color = "black",
    stroke = 1.1
  ) +
  geom_text(
    data = unavailable,
    label = "/",
    size = 3.2,
    color = "black",
    fontface = "bold"
  ) +
  geom_point(
    data = insufficient,
    shape = 3,
    size = 2,
    color = "black",
    stroke = 1.1
  ) +
  facet_grid(. ~ Dataset, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = spectral, limits = c(-1, 1)) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Pearson",
    caption = "x: constant prediction; /: unavailable or no reference; +: insufficient samples"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.45),
    strip.background = element_rect(fill = "#F7E8EA", color = "black", linewidth = 0.45),
    strip.text = element_text(size = 12, color = "#7F2725", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.caption = element_text(hjust = 0, color = "black"),
    legend.position = "right"
  )

ggsave(output_file, p, width = 11.69, height = 8.27, units = "in")
