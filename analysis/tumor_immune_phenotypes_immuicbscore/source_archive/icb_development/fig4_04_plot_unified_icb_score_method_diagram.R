suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

if (!requireNamespace("png", quietly = TRUE)) {
  stop("Package 'png' is required. Please install it with install.packages('png').")
}

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
out_dir <- file.path(fig4_dir, "Unified_ICB_score_cluster_ICB_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

material_file <- file.path(out_dir, "Unified_ICB_score_method_GPT_image_material_less_pies.png")
fallback_file <- file.path(out_dir, "Unified_ICB_score_method_GPT_image_material_detailed_RF.png")
if (!file.exists(material_file) && file.exists(fallback_file)) {
  material_file <- fallback_file
}
stopifnot(file.exists(material_file))
img <- png::readPNG(material_file)

theme_method <- theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 12, 8, 12)
  )

step_df <- data.frame(
  x = c(1.20, 3.20, 5.20, 7.30, 9.25),
  y = rep(5.08, 5),
  w = c(1.55, 1.52, 1.70, 1.95, 1.55),
  h = rep(0.58, 5),
  fill = c("#F7E7EA", "#E8F1F2", "#E8F1F2", "#FFF4D6", "#EAF4EA"),
  label = c(
    "ImmuCellAI2\n53 fractions",
    "CLR\ntransform",
    "Training-fold\nscaling",
    "Balanced ranger\nrandom forest",
    "Unified ICB\nScore"
  ),
  stringsAsFactors = FALSE
)
step_df$xmin <- step_df$x - step_df$w / 2
step_df$xmax <- step_df$x + step_df$w / 2
step_df$ymin <- step_df$y - step_df$h / 2
step_df$ymax <- step_df$y + step_df$h / 2

detail_df <- data.frame(
  x = c(1.20, 3.20, 5.20, 7.30, 9.25),
  y = c(3.22, 3.22, 3.22, 3.18, 3.22),
  w = c(1.70, 1.80, 1.90, 2.45, 1.78),
  h = c(0.98, 1.14, 1.14, 1.38, 0.98),
  label = c(
    "p_i = (p_i1,...,p_i53)\nimmune-cell fractions only",
    "CLR_ij = log(p_ij + eps)\n - mean_j log(p_ij + eps)\neps = 1e-5",
    "z_ij = (CLR_ij - mu_j,train)\n / sigma_j,train\ntraining fold only",
    "class weights:\nw_1 = 0.5 / Pr(y=1)\nw_0 = 0.5 / Pr(y=0)\n1200 trees; mtry = 7\nmin.node.size = 8",
    "Score_i = Pr(y_i = 1 | x_i)\nhigher = more likely responder"
  ),
  stringsAsFactors = FALSE
)
detail_df$xmin <- detail_df$x - detail_df$w / 2
detail_df$xmax <- detail_df$x + detail_df$w / 2
detail_df$ymin <- detail_df$y - detail_df$h / 2
detail_df$ymax <- detail_df$y + detail_df$h / 2

arrow_df <- data.frame(
  x = c(1.95, 3.95, 6.02, 8.20),
  y = rep(5.08, 4),
  xend = c(2.40, 4.34, 6.46, 8.55),
  yend = rep(5.08, 4)
)

rf_note <- data.frame(
  x = c(6.62, 7.30, 7.98),
  y = c(2.02, 2.02, 2.02),
  label = c("bootstrap\nsamples", "random feature\nsplits", "probability\nvote")
)

p <- ggplot() +
  annotation_raster(img, xmin = 0.45, xmax = 10.05, ymin = 1.00, ymax = 4.58, interpolate = TRUE) +
  geom_rect(
    data = detail_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    color = "#2C2C2C",
    linewidth = 0.32,
    alpha = 0.88
  ) +
  geom_rect(
    data = step_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    color = "black",
    linewidth = 0.43,
    alpha = 0.98
  ) +
  scale_fill_identity() +
  geom_segment(
    data = arrow_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.13, "inches"), type = "closed"),
    linewidth = 0.52,
    lineend = "round",
    color = "grey18"
  ) +
  geom_text(
    data = step_df,
    aes(x = x, y = y, label = label),
    size = 3.04,
    fontface = "bold",
    lineheight = 0.90,
    color = "black"
  ) +
  geom_text(
    data = detail_df,
    aes(x = x, y = y, label = label),
    size = 2.08,
    color = "grey13",
    lineheight = 0.86
  ) +
  geom_point(
    data = rf_note,
    aes(x = x, y = y),
    shape = 21,
    size = 2.70,
    stroke = 0.42,
    color = "black",
    fill = "#FAC74C"
  ) +
  geom_text(
    data = rf_note,
    aes(x = x, y = y - 0.34, label = label),
    size = 1.90,
    color = "grey18",
    lineheight = 0.82
  ) +
  annotate(
    "text",
    x = 5.20,
    y = 5.90,
    label = "Unified ICB Score calculation",
    size = 4.50,
    fontface = "bold",
    color = "black"
  ) +
  annotate(
    "text",
    x = 5.20,
    y = 5.55,
    label = "CLR-transformed 53 ImmuCellAI2 fractions + balanced ranger random forest",
    size = 2.92,
    color = "grey25"
  ) +
  annotate(
    "text",
    x = 5.20,
    y = 0.35,
    label = "Only immune-cell fractions are used as model inputs; downstream validation modules are omitted from this method-only diagram.",
    size = 2.35,
    color = "grey25"
  ) +
  coord_cartesian(xlim = c(0.25, 10.15), ylim = c(0.12, 6.10), expand = FALSE) +
  theme_method

ggsave(
  file.path(out_dir, "Figure0_Unified_ICB_score_method_diagram.pdf"),
  p, width = 10.2, height = 5.8, units = "in"
)
ggsave(
  file.path(out_dir, "Figure0_Unified_ICB_score_method_diagram.png"),
  p, width = 10.2, height = 5.8, units = "in", dpi = 420
)

method_text <- c(
  "Unified ICB Score calculation, GPT-image style with fewer pie-chart elements",
  "",
  "The method diagram uses a no-text GPT-image visual material with sparse decorative pie glyphs.",
  "All scientific labels, formulas, and model parameters are overlaid by R.",
  "",
  "Step 1. ImmuCellAI2 estimates a 53-cell fraction vector:",
  "  p_i = (p_i1,...,p_i53).",
  "Step 2. Fractions are CLR-transformed:",
  "  CLR_ij = log(p_ij + eps) - mean_j[log(p_ij + eps)], eps = 1e-5.",
  "Step 3. In each training fold, features are standardized:",
  "  z_ij = (CLR_ij - mu_j,train) / sigma_j,train.",
  "Step 4. Balanced ranger random forest:",
  "  w_1 = 0.5 / Pr(y=1), w_0 = 0.5 / Pr(y=0);",
  "  num.trees = 1200; mtry = 7; min.node.size = 8.",
  "Step 5. Unified ICB Score:",
  "  Score_i = Pr(y_i = 1 | x_i)."
)
writeLines(method_text, file.path(out_dir, "Unified_ICB_score_method_description.txt"), useBytes = TRUE)

message("Unified ICB Score method diagram with fewer pie elements written to: ", out_dir)
