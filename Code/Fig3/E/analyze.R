# Preserve measured times; derive the log-scale plotting coordinate in R.
d <- read.delim(file.path(data_dir, "E", "runtime_plot_data.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
stopifnot(nrow(d) == 14L, !anyDuplicated(d[c("MethodShort", "Scenario")]),
          setequal(d$MethodShort, methods), all(d$Status == "OK"),
          setequal(d$Scenario, c("single_sample_1core", "fifty_samples_8cores")),
          all(d$Samples == ifelse(d$Scenario == "single_sample_1core", 1L, 50L)),
          all(d$Cores == ifelse(d$Scenario == "single_sample_1core", 1L, 8L)),
          all(is.finite(d$SecondsPerSample)), all(d$SecondsPerSample >= 0))
# ElapsedSeconds is stored at lower precision; allow its rounding error.
same_numbers(d$ElapsedSeconds / d$Samples, d$SecondsPerSample, tolerance = 1e-7)
log_time <- log10(d$SecondsPerSample + 1)
same_numbers(log_time, d$LogTime)
d$LogTime <- log_time
analysis$E <- d
