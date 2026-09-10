# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

work_dir <- fig3_work_dir()
fig2_dir <- file.path(work_dir, "Fig2")
out_dir <- file.path(fig2_dir, "mean_pearson_compare")
pbmc_dir <- file.path(fig2_dir, "Newman_Monaco_7tools_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

four_dataset_file <- file.path(out_dir, "overall_summary_current_plot_scope.txt")
newman_file <- file.path(
  pbmc_dir, "RNA_Sieve_newman_pbmcs", "seven_tools_compare",
  "comparison_summary_7tools.txt"
)
monaco_file <- file.path(
  pbmc_dir, "RNA_Sieve_monaco_pbmcs", "seven_tools_compare",
  "comparison_summary_7tools.txt"
)

required_files <- c(four_dataset_file, newman_file, monaco_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required.")
library(ggplot2)

method_map <- c(
  "ImmuCellAI 2.0" = "ImmuCellAI 2.0",
  ImmuCellAI2_flat_VB = "ImmuCellAI 2.0",
  ImmuCellAI2_tcell_VB_UNKNOWN = "ImmuCellAI 2.0",
  BayesPrism = "BayesPrism",
  BayesPrism_first_state_chain600_burn500 = "BayesPrism",
  CIBERSORT = "CIBERSORT",
  CIBERSORT_default_LM22 = "CIBERSORT",
  CITMIC = "CITMIC",
  CITMIC_native = "CITMIC",
  DWLS = "DWLS",
  DWLS_native_ref_top1000 = "DWLS",
  ImmuCellAI = "ImmuCellAI",
  ImmuCellAI_native = "ImmuCellAI",
  MuSiC = "MuSiC",
  MuSiC_native_music.basic = "MuSiC"
)

method_order <- c(
  "ImmuCellAI 2.0", "BayesPrism", "CIBERSORT", "CITMIC",
  "DWLS", "ImmuCellAI", "MuSiC"
)
dataset_order <- c(
  "GSE146771", "GSE164522", "GSE176078",
  "GSE107011", "Newman PBMC", "Monaco PBMC"
)

four <- read.delim(four_dataset_file, check.names = FALSE)
four <- four[, c("Dataset", "Method", "N_targets", "MeanPearson", "ValidPearsonTargets"), drop = FALSE]

read_pbmc <- function(file, dataset) {
  x <- read.delim(file, check.names = FALSE)
  x$Dataset <- dataset
  x[, c("Dataset", "Method", "N_targets", "MeanPearson", "ValidPearsonTargets"), drop = FALSE]
}

pbmc <- rbind(
  read_pbmc(newman_file, "Newman PBMC"),
  read_pbmc(monaco_file, "Monaco PBMC")
)
plot_data <- rbind(four, pbmc)
plot_data$MethodOriginal <- as.character(plot_data$Method)
plot_data$Method <- unname(method_map[plot_data$MethodOriginal])
plot_data$MeanPearson <- suppressWarnings(as.numeric(plot_data$MeanPearson))
plot_data$N_targets <- suppressWarnings(as.numeric(plot_data$N_targets))
plot_data$ValidPearsonTargets <- suppressWarnings(as.numeric(plot_data$ValidPearsonTargets))
plot_data <- plot_data[
  plot_data$Method %in% method_order & plot_data$Dataset %in% dataset_order,
  , drop = FALSE
]
plot_data$Dataset <- factor(plot_data$Dataset, levels = dataset_order)
plot_data$Method <- factor(plot_data$Method, levels = method_order)
plot_data <- plot_data[order(plot_data$Dataset, plot_data$Method), , drop = FALSE]

write.table(
  plot_data,
  file.path(out_dir, "six_datasets_mean_pearson_plot_data.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

method_colors <- c(
  "ImmuCellAI 2.0" = "#F66463",
  BayesPrism = "#379DA5",
  CIBERSORT = "#AE997E",
  CITMIC = "#95A8AC",
  DWLS = "#6BA5C3",
  ImmuCellAI = "#8E6BBE",
  MuSiC = "#FAC74C"
)

plot_theme <- theme_bw(base_size = 9) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "gray85", linewidth = 0.35, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.20),
    strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 0.20),
    strip.text = element_text(size = 8.2, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.3, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.12),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 6.5, colour = "black"),
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(5, "mm"),
    plot.margin = margin(3, 4, 3, 3, "mm")
  )

p <- ggplot(plot_data, aes(x = Method, y = MeanPearson, fill = Method)) +
  geom_col(width = 0.70, color = "#202020", linewidth = 0.16) +
  facet_wrap(~ Dataset, ncol = 6) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.035))
  ) +
  labs(x = NULL, y = "Mean Pearson correlation") +
  plot_theme

ggsave(
  file.path(out_dir, "overall_mean_pearson_bar.pdf"),
  p, width = 11.69, height = 3.35, units = "in", bg = "white", limitsize = FALSE
)
ggsave(
  file.path(out_dir, "overall_mean_pearson_bar.png"),
  p, width = 11.69, height = 3.35, units = "in", dpi = 500, bg = "white", limitsize = FALSE
)

print(plot_data[, c("Dataset", "Method", "N_targets", "ValidPearsonTargets", "MeanPearson")])
message("Saved six-dataset mean Pearson plot to: ", out_dir)
