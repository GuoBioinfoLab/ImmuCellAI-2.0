# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

workspace_dir <- fig3_work_dir()
out_dir <- file.path(workspace_dir, "Fig2", "runtime_compare")
summary_file <- file.path(out_dir, "runtime_benchmark_1sample_vs_50samples_summary.txt")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required.")
}

library(ggplot2)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(summary_file) &&
    Sys.getenv("IMMUCELLAI2_FIG3_ARCHIVED_TIMING", "0") != "1") {
  stop("Runtime summary missing. Run the timing benchmark or explicitly enable archived timing constants.")
}
if (file.exists(summary_file)) {
  combined <- read.delim(summary_file, check.names = FALSE)
} else {
  combined <- data.frame(
    Method = rep(c(
      "BayesPrism_first_state_chain600_burn500",
      "CIBERSORT_nuSVR_ref_top1000",
      "CITMIC_native",
      "DWLS_weighted_lm",
      "ImmuCellAI_native",
      "ImmuCellAI2_tcell_VB",
      "MuSiC_basic"
    ), each = 2),
    Status = "OK",
    ElapsedSeconds = c(
      161.6405311, 2492.3269670,
      1.5447578, 129.4113412,
      4.0495570, 67.7533329,
      1.3967259, 6.6229210,
      0.7796760, 1.8751221,
      0.3477101, 4.9696131,
      8.4314771, 40.0910151
    ),
    SecondsPerSample = c(
      161.64053106, 49.84653934,
      1.54475784, 2.58822682,
      4.04955697, 1.35506666,
      1.39672589, 0.13245842,
      0.77967596, 0.03750244,
      0.34771013, 0.09939226,
      8.43147707, 0.80182030
    ),
    Samples = rep(c(1, 50), 7),
    GenesUsed = rep(c(59073, 1000, 59073, 59073, 59073, 5505, 59073), each = 2),
    CellStatesReturned = rep(c(53, 53, 86, 53, 25, 53, 53), each = 2),
    Cores = rep(c(1, 8), 7),
    Scenario = rep(c("single_sample_1core", "fifty_samples_8cores"), 7),
    stringsAsFactors = FALSE
  )
  write.table(
    combined,
    summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
combined <- combined[combined$Status == "OK", , drop = FALSE]
combined$Scenario <- factor(
  combined$Scenario,
  levels = c("single_sample_1core", "fifty_samples_8cores"),
  labels = c("1 sample / 1 core", "50 samples / 8 cores")
)

combined$Log10ElapsedPlus1 <- log10(combined$ElapsedSeconds + 1)
combined$Log10SecondsPerSamplePlus1 <- log10(combined$SecondsPerSample + 1)

method_order_top_to_bottom <- c(
  "ImmuCellAI2_tcell_VB",
  "ImmuCellAI_native",
  "DWLS_weighted_lm",
  "CIBERSORT_nuSVR_ref_top1000",
  "CITMIC_native",
  "MuSiC_basic",
  "BayesPrism_first_state_chain600_burn500"
)
method_order_top_to_bottom <- method_order_top_to_bottom[
  method_order_top_to_bottom %in% unique(combined$Method)
]
combined$Method <- factor(combined$Method, levels = rev(method_order_top_to_bottom))

my_cols <- c("1 sample / 1 core" = "#379DA5", "50 samples / 8 cores" = "#F66463")

runtime_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.7),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    legend.title = element_blank()
  )

p_elapsed <- ggplot(combined, aes(x = Method, y = Log10ElapsedPlus1, fill = Scenario)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66,
           color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = NULL,
    y = "log10(elapsed seconds + 1)",
    title = "Runtime comparison"
  ) +
  runtime_theme

p_per_sample <- ggplot(combined, aes(x = Method, y = Log10SecondsPerSamplePlus1, fill = Scenario)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66,
           color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = NULL,
    y = "log10(seconds per sample + 1)",
    title = "Per-sample runtime comparison"
  ) +
  runtime_theme

write.table(
  combined,
  file.path(out_dir, "runtime_benchmark_1sample_vs_50samples_logvalue_plot_data.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(file.path(out_dir, "runtime_elapsed_1sample_vs_50samples_logvalue.png"),
       p_elapsed, width = 11, height = 5.6, dpi = 300)
ggsave(file.path(out_dir, "runtime_elapsed_1sample_vs_50samples_logvalue.pdf"),
       p_elapsed, width = 11, height = 5.6)
ggsave(file.path(out_dir, "runtime_seconds_per_sample_1sample_vs_50samples_logvalue.png"),
       p_per_sample, width = 11, height = 5.6, dpi = 300)
ggsave(file.path(out_dir, "runtime_seconds_per_sample_1sample_vs_50samples_logvalue.pdf"),
       p_per_sample, width = 11, height = 5.6)

print(combined[, c("Method", "Scenario", "ElapsedSeconds", "SecondsPerSample",
                   "Log10ElapsedPlus1", "Log10SecondsPerSamplePlus1")])
