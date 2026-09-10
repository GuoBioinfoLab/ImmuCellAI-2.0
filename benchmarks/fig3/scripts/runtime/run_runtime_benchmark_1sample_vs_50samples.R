# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

workspace_dir <- fig3_work_dir()
runtime_compare_dir <- file.path(workspace_dir, "Fig2", "runtime_compare")
benchmark_script <- fig3_script("runtime/benchmark_7tools_runtime_efficiency.R")
out_dir <- runtime_compare_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stable_tmp <- file.path(workspace_dir, "Rtmp_immucellai2")
dir.create(stable_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = stable_tmp, TEMP = stable_tmp, TMPDIR = stable_tmp)

if (!file.exists(benchmark_script)) {
  stop("Cannot find benchmark script: ", benchmark_script)
}

run_scenario <- function(label, samples, cores) {
  message("Running runtime scenario: ", label)
  env <- c(
    IMMUCELLAI2_RUNTIME_SAMPLES = as.character(samples),
    IMMUCELLAI2_RUNTIME_CORES = as.character(cores),
    IMMUCELLAI2_RUNTIME_LABEL = label,
    TMP = stable_tmp,
    TEMP = stable_tmp,
    TMPDIR = stable_tmp
  )
  old_env <- Sys.getenv(names(env), unset = NA_character_)
  do.call(Sys.setenv, as.list(env))
  on.exit({
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, setNames(list(old_env[[nm]]), nm))
      }
    }
  }, add = TRUE)
  log_file <- file.path(out_dir, paste0(label, "_Rscript_stdout.log"))
  err_file <- file.path(out_dir, paste0(label, "_Rscript_stderr.log"))
  status <- system2(
    file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
    args = shQuote(benchmark_script),
    stdout = log_file,
    stderr = err_file
  )
  if (!identical(status, 0L)) {
    stop("Scenario failed: ", label, ". See logs: ", log_file, " and ", err_file)
  }
  summary_file <- file.path(runtime_compare_dir, paste0("runtime_benchmark_7tools_", label),
                            "runtime_benchmark_7tools_summary.txt")
  if (!file.exists(summary_file)) {
    stop("Scenario summary not found: ", summary_file)
  }
  x <- read.delim(summary_file, check.names = FALSE)
  x$Scenario <- label
  x$ScenarioSamples <- samples
  x$ScenarioCores <- cores
  x
}

scenario_1 <- run_scenario("single_sample_1core", samples = 1, cores = 1)
scenario_2 <- run_scenario("fifty_samples_8cores", samples = 50, cores = 8)

combined <- rbind(scenario_1, scenario_2)
combined$ThroughputSamplesPerSecond <- combined$Samples / combined$ElapsedSeconds
combined$Scenario <- factor(combined$Scenario,
                            levels = c("single_sample_1core", "fifty_samples_8cores"))
combined <- combined[order(combined$Method, combined$Scenario), ]

write.table(
  combined,
  file.path(out_dir, "runtime_benchmark_1sample_vs_50samples_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

wide <- reshape(
  combined[, c("Method", "Scenario", "ElapsedSeconds", "SecondsPerSample",
               "ThroughputSamplesPerSecond", "Status")],
  idvar = "Method",
  timevar = "Scenario",
  direction = "wide"
)
write.table(
  wide,
  file.path(out_dir, "runtime_benchmark_1sample_vs_50samples_wide.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  ok <- combined[combined$Status == "OK", , drop = FALSE]
  ok$Method <- factor(ok$Method, levels = rev(unique(ok$Method[order(ok$ElapsedSeconds)])))

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
