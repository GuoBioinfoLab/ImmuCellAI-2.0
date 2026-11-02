# Sourced by ../reproduce_Fig4.R in panel order.
# A: keep the exact saved row/column order, not a new randomized heatmap.
z <- read_tab("A/heatmap_plot_data.tsv")
rownames(z) <- z$cell_type
z$cell_type <- NULL
z <- as.matrix(z)
raw <- read_tab("A/cell_fractions.tsv")
rownames(raw) <- raw$sample
raw$sample <- NULL
ann <- read_tab("A/sample_annotation.tsv")
check("A_53_by_4854", identical(dim(z), c(53L,4854L)))
check("A_order_unique", !anyDuplicated(ann$sample) && identical(ann$sample, colnames(z)))
z_test <- t(scale(as.matrix(raw[colnames(z), rownames(z)])))
check("A_z_scores_reproduced", near(z, pmax(-1, pmin(1, z_test))))
