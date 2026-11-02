if (!is.null(run_dir)) stop("Runtime is a separate historical experiment.")
d <- pipeline_read(cfg$timings)
stopifnot(nrow(d) == 14L, all(d$Status == "OK"), all(d$Samples > 0),
          all(is.finite(d$ElapsedSeconds)), all(d$ElapsedSeconds >= 0),
          !anyDuplicated(d[c("MethodShort", "Scenario")]),
          all(abs(d$SecondsPerSample - d$ElapsedSeconds / d$Samples) < 1e-7))
d$LogTime <- log10(d$SecondsPerSample + 1)
d
