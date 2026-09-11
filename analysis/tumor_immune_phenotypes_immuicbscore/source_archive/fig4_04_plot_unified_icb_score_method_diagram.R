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

material_file <- file.path(out_dir, "Unified_ICB_score_method_GPT_image_material.png")
stopifnot(file.exists(material_file))

img <- png::readPNG(material_file)

theme_method <- theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 12, 10, 12)
  )

label_box <- data.frame(
  x = c(1.25, 3.30, 5.35, 7.40, 9.25),
  y = c(5.18, 5.18, 5.18, 5.18, 5.18),
  w = c(1.55, 1.65, 1.65, 1.75, 1.55),
  h = c(0.56, 0.56, 0.56, 0.56, 0.56),
  fill = c("#F7E7EA", "#E8F1F2", "#E8F1F2", "#E8F1F2", "#EAF4EA"),
  label = c(
    "53-cell fractions",
    "CLR transform",
    "Training-fold\nstandardization",
    "Balanced ranger\nrandom forest",
    "Unified ICB\nScore"
  ),
  stringsAsFactors = FALSE
)
label_box$xmin <- label_box$x - label_box$w / 2
label_box$xmax <- label_box$x + label_box$w / 2
label_box$ymin <- label_box$y - label_box$h / 2
label_box$ymax <- label_box$y + label_box$h / 2

formula_df <- data.frame(
  x = c(3.30, 5.35, 7.40, 9.25),
  y = c(4.42, 4.42, 4.42, 4.42),
  label = c(
    "CLR_ij = log(p_ij + 1e-5)\n          - mean_j log(p_ij + 1e-5)",
    "z_ij = (CLR_ij - mu_j,train)\n       / sigma_j,train",
    "num.trees = 1200\nmtry = floor(sqrt(53))\nmin.node.size = 8\nclass weight = 0.5 / prevalence",
    "Score_i = Pr(response_i = 1)"
  ),
  stringsAsFactors = FALSE
)

validation_box <- data.frame(
  x = c(8.05, 9.8),
  y = c(1.12, 1.12),
  w = c(1.65, 1.55),
  h = c(0.56, 0.56),
  fill = c("#FFF4D6", "#FFF4D6"),
  label = c("TCGA C1-C4\nassociation", "External ICB\nROC validation"),
  stringsAsFactors = FALSE
)
validation_box$xmin <- validation_box$x - validation_box$w / 2
validation_box$xmax <- validation_box$x + validation_box$w / 2
validation_box$ymin <- validation_box$y - validation_box$h / 2
validation_box$ymax <- validation_box$y + validation_box$h / 2

arrow_df <- data.frame(
  x = c(1.96, 4.08, 6.12, 8.20, 8.93, 9.62),
  y = c(5.18, 5.18, 5.18, 5.18, 4.82, 4.82),
  xend = c(2.46, 4.58, 6.67, 8.62, 8.25, 9.68),
  yend = c(5.18, 5.18, 5.18, 5.18, 1.50, 1.50),
  stringsAsFactors = FALSE
)

caption_df <- data.frame(
  x = 5.2,
  y = 0.42,
  label = "Model input is restricted to ImmuCellAI2-derived immune-cell fractions; TCGA cluster and ICB ROC analyses are downstream validations."
)

p <- ggplot() +
  annotation_raster(img, xmin = 0.55, xmax = 10.35, ymin = 1.55, ymax = 4.15, interpolate = TRUE) +
  geom_rect(
    data = label_box,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    color = "black",
    linewidth = 0.42,
    alpha = 0.96
  ) +
  geom_rect(
    data = validation_box,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    color = "black",
    linewidth = 0.42,
    alpha = 0.96
  ) +
  scale_fill_identity() +
  geom_segment(
    data = arrow_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.14, "inches"), type = "closed"),
    linewidth = 0.5,
    lineend = "round",
    color = "grey20"
  ) +
  geom_text(
    data = label_box,
    aes(x = x, y = y, label = label),
    size = 3.0,
    fontface = "bold",
    lineheight = 0.92,
    color = "black"
  ) +
  geom_text(
    data = validation_box,
    aes(x = x, y = y, label = label),
    size = 2.9,
    fontface = "bold",
    lineheight = 0.92,
    color = "black"
  ) +
  geom_text(
    data = formula_df,
    aes(x = x, y = y, label = label),
    size = 2.35,
    color = "grey20",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 5.2,
    y = 5.92,
    label = "Unified ICB Score calculation",
    size = 4.35,
    fontface = "bold",
    color = "black"
  ) +
  geom_text(
    data = caption_df,
    aes(x = x, y = y, label = label),
    size = 2.45,
    color = "grey25"
  ) +
  coord_cartesian(xlim = c(0.25, 10.6), ylim = c(0.20, 6.15), expand = FALSE) +
  theme_method

ggsave(
  file.path(out_dir, "Figure0_Unified_ICB_score_method_diagram.pdf"),
  p, width = 10.4, height = 6.0, units = "in"
)
ggsave(
  file.path(out_dir, "Figure0_Unified_ICB_score_method_diagram.png"),
  p, width = 10.4, height = 6.0, units = "in", dpi = 360
)

method_text <- c(
  "Unified ICB Score calculation with GPT-image material figure",
  "",
  "The final method figure uses a GPT-image generated no-text biomedical illustration as visual material.",
  "All scientific labels, equations, model parameters, and validation labels are overlaid by R to avoid image-generated text errors.",
  "",
  "Step 1. ImmuCellAI2 estimates a 53-cell fraction vector for each sample.",
  "Step 2. Fractions are transformed by centered log-ratio transformation:",
  "  CLR_ij = log(p_ij + 1e-5) - mean_j[log(p_ij + 1e-5)].",
  "Step 3. In each training fold, features are standardized as:",
  "  z_ij = (CLR_ij - mu_j,train) / sigma_j,train.",
  "Step 4. A balanced ranger probability random forest is trained:",
  "  num.trees = 1200;",
  "  mtry = floor(sqrt(53));",
  "  min.node.size = 8;",
  "  class weight = 0.5 / class prevalence.",
  "Step 5. Unified ICB Score is the predicted responder probability:",
  "  Score_i = Pr(response_i = 1).",
  "",
  "TCGA C1-C4 association and external ICB ROC analyses are downstream validation steps, not model inputs."
)
writeLines(method_text, file.path(out_dir, "Unified_ICB_score_method_description.txt"), useBytes = TRUE)

message("GPT-image enhanced method diagram written to: ", out_dir)
