# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

workspace_dir <- fig3_work_dir()
out_dir <- file.path(workspace_dir, "Fig2", "runtime_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario_info <- data.frame(
  Scenario = c("single_sample_1core", "fifty_samples_8cores"),
  ScenarioSamples = c(1, 50),
  ScenarioCores = c(1, 8),
  stringsAsFactors = FALSE
)

read_scenario <- function(i) {
  label <- scenario_info$Scenario[i]
  summary_file <- file.path(out_dir, paste0("runtime_benchmark_7tools_", label),
                            "runtime_benchmark_7tools_summary.txt")
  if (!file.exists(summary_file)) {
    stop("Missing scenario summary: ", summary_file)
  }
  x <- read.delim(summary_file, check.names = FALSE)
  x$Scenario <- label
  x$ScenarioSamples <- scenario_info$ScenarioSamples[i]
  x$ScenarioCores <- scenario_info$ScenarioCores[i]
  x
}

combined <- do.call(rbind, lapply(seq_len(nrow(scenario_info)), read_scenario))
combined$ThroughputSamplesPerSecond <- combined$Samples / combined$ElapsedSeconds
combined$Scenario <- factor(combined$Scenario,
                            levels = c("single_sample_1core", "fifty_samples_8cores"))
combined <- combined[order(combined$Method, combined$Scenario), ]

write.table(combined,
            file.path(out_dir, "runtime_benchmark_1sample_vs_50samples_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

wide <- reshape(
  combined[, c("Method", "Scenario", "Status", "ElapsedSeconds", "SecondsPerSample",
               "ThroughputSamplesPerSecond", "Samples", "Cores", "GenesUsed",
               "CellStatesReturned")],
  idvar = "Method",
  timevar = "Scenario",
  direction = "wide"
)
write.table(wide,
            file.path(out_dir, "runtime_benchmark_1sample_vs_50samples_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  ok <- combined[combined$Status == "OK", , drop = FALSE]
  method_order <- unique(ok$Method[order(ok$SecondsPerSample)])
  ok$Method <- factor(ok$Method, levels = rev(method_order))

  p_elapsed <- ggplot(ok, aes(x = Method, y = ElapsedSeconds, fill = Scenario)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.66,
             color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_y_log10() +
    scale_fill_manual(values = c(single_sample_1core = "#379DA5",
                                 fifty_samples_8cores = "#F66463")) +
    labs(x = NULL, y = "Elapsed seconds, log10 scale",
         title = "Runtime comparison: 1 sample / 1 core vs 50 samples / 8 cores") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 0.7),
          plot.title = element_text(face = "bold"))

  p_per_sample <- ggplot(ok, aes(x = Method, y = SecondsPerSample, fill = Scenario)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.66,
             color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_y_log10() +
    scale_fill_manual(values = c(single_sample_1core = "#379DA5",
                                 fifty_samples_8cores = "#F66463")) +
    labs(x = NULL, y = "Seconds per sample, log10 scale",
         title = "Per-sample runtime comparison") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 0.7),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(out_dir, "runtime_elapsed_1sample_vs_50samples.png"),
         p_elapsed, width = 11, height = 5.6, dpi = 300)
  ggsave(file.path(out_dir, "runtime_elapsed_1sample_vs_50samples.pdf"),
         p_elapsed, width = 11, height = 5.6)
  ggsave(file.path(out_dir, "runtime_seconds_per_sample_1sample_vs_50samples.png"),
         p_per_sample, width = 11, height = 5.6, dpi = 300)
  ggsave(file.path(out_dir, "runtime_seconds_per_sample_1sample_vs_50samples.pdf"),
         p_per_sample, width = 11, height = 5.6)
}

print(combined)
