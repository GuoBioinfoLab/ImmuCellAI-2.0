read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             fileEncoding = "UTF-8", quote = "", comment.char = "")
}

key <- function(dat, fields) do.call(paste, c(dat[fields], sep = "|"))

verify_data <- function(data_dir) {
  largest <- 0
  comparisons <- 0L
  equal_num <- function(x, y, label) {
    stopifnot(length(x) == length(y), identical(is.na(x), is.na(y)))
    difference <- if (all(is.na(x))) 0 else max(abs(x - y), na.rm = TRUE)
    if (!is.finite(difference) || difference > 1e-12) stop(label, ": difference=", difference)
    largest <<- max(largest, difference)
    comparisons <<- comparisons + 1L
  }
  same_table <- function(actual, expected, keys, fields, label) {
    ka <- key(actual, keys)
    ke <- key(expected, keys)
    stopifnot(!anyDuplicated(ka), !anyDuplicated(ke), setequal(ka, ke))
    expected <- expected[match(ka, ke), , drop = FALSE]
    for (field in fields) {
      if (is.numeric(actual[[field]])) equal_num(actual[[field]], expected[[field]], paste(label, field))
      else if (is.logical(expected[[field]])) stopifnot(identical(as.logical(actual[[field]]), expected[[field]]))
      else if (!identical(as.character(actual[[field]]), as.character(expected[[field]]))) stop(label, ": text mismatch in ", field)
    }
  }
  get <- function(path) read_tsv(file.path(data_dir, path))
  clinical <- get("_shared/clinical_sheet1.tsv")
  workbook <- as.data.frame(readxl::read_excel(file.path(data_dir, "_shared/original_metadata/PRJNA683803.xlsx")))
  same_table(clinical, workbook, "Run", setdiff(names(clinical), "Run"), "clinical workbook")
  original_meta <- read.csv(file.path(data_dir, "_shared/original_metadata/PRJNA683803_combined.csv"), check.names = FALSE)
  stopifnot(!anyDuplicated(original_meta$Run), setequal(original_meta$Run, clinical$Run))
  original_meta <- original_meta[match(clinical$Run, original_meta$Run), ]
  stopifnot(identical(tolower(original_meta$outcome), tolower(clinical$outcome)))
  mapping <- get("_shared/time_mapping.tsv")
  stopifnot(!anyDuplicated(mapping$grouping))
  clinical$Time <- mapping$Time[match(clinical$grouping, mapping$grouping)]
  clinical$TimeDay <- mapping$TimeDay[match(clinical$grouping, mapping$grouping)]
  outcomes <- c(died = "Died", dead = "Died", survived = "Survived", survive = "Survived",
                "died-iris" = "Died-IRIS", "dead-iris" = "Died-IRIS",
                "survived-iris" = "Survived-IRIS", "survive-iris" = "Survived-IRIS")
  clinical$status <- unname(outcomes[tolower(trimws(clinical$outcome))])
  clinical$included <- clinical$status %in% c("Died", "Survived") & !is.na(clinical$Time)
  same_table(get("_shared/sample_eligibility.tsv"), clinical, "Run",
             c("Patient_ID", "status", "Time", "TimeDay", "included"), "eligibility")
  eligible <- clinical[clinical$included, ]
  stopifnot(nrow(clinical) == 211, nrow(eligible) == 184, !anyDuplicated(clinical$Run),
            length(unique(eligible$Patient_ID)) == 183)
  fractions <- get("_shared/state_fraction.tsv")
  names(fractions)[1] <- "Run"
  stopifnot(nrow(fractions) == 211, ncol(fractions) == 54, !anyDuplicated(fractions$Run),
            setequal(fractions$Run, clinical$Run), all(is.finite(as.matrix(fractions[-1]))),
            all(as.matrix(fractions[-1]) >= 0), all(as.matrix(fractions[-1]) <= 1))
  equal_num(rowSums(fractions[-1]), rep(1, 211), "53-state proportions sum")
  combined <- get("_shared/state_fraction_with_metadata.tsv")
  same_table(combined, fractions, "Run", names(fractions)[-1], "combined fraction table")
  index <- get("panel_index.tsv")
  stopifnot(identical(index$panel, c("A", "B")))
  panels <- list()
  for (panel in c("A", "B")) {
    cells <- get(paste0(panel, "/cell_order.tsv"))$CellType
    expected_cells <- if (panel == "A") c("CD4Tn", "CD8Tn", "Th1", "Tr1") else c("Tc", "cNK", "Neutrophil", "MDSC", "Th1", "CD4Tn", "CD8Tn", "Tr1")
    stopifnot(identical(cells, expected_cells))
    samples <- get(paste0(panel, "/sample_level.tsv"))
    summary <- get(paste0(panel, "/plot_data.tsv"))
    stopifnot(nrow(samples) == 184 * length(cells), nrow(summary) == 6 * length(cells),
              !anyDuplicated(key(samples, c("Run", "CellType"))), setequal(samples$Run, eligible$Run),
              identical(unique(samples$CellType), cells))
    matched <- eligible[match(samples$Run, eligible$Run), ]
    stopifnot(identical(samples$status, matched$status), identical(samples$Time, matched$Time))
    if (panel == "B") {
      stopifnot(identical(samples$Patient_ID, matched$Patient_ID))
      equal_num(samples$TimeDay, matched$TimeDay, "numeric time")
    }
    values <- vapply(seq_len(nrow(samples)), function(i) fractions[[samples$CellType[i]]][match(samples$Run[i], fractions$Run)], numeric(1))
    equal_num(samples$Abundance, values, paste(panel, "frozen sample values"))
    original_samples <- get(paste0("_shared/", panel, "_original_samples.tsv"))
    same_table(samples, original_samples, c("Run", "CellType"), c("status", "Time", "Abundance"), paste(panel, "source samples"))
    original_summary <- get(if (panel == "A") "_shared/A_original_means.tsv" else "_shared/B_original_summary.tsv")
    same_table(summary, original_summary, c("CellType", "status", "Time"), if (panel == "A") "Mean" else c("Mean", "SEM", "N", "TimeDay"), paste(panel, "source summary"))
    for (i in seq_len(nrow(summary))) {
      row <- summary[i, ]
      x <- samples$Abundance[samples$CellType == row$CellType & samples$status == row$status & samples$Time == row$Time]
      equal_num(row$Mean, mean(x), paste(panel, row$CellType, "mean"))
      if (panel == "B") equal_num(c(row$N, row$SEM), c(length(x), sd(x) / sqrt(length(x))), paste(panel, row$CellType, "SEM"))
    }
    counts <- get(paste0(panel, "/sample_counts.tsv"))
    stopifnot(nrow(counts) == length(cells) * 8, !anyDuplicated(key(counts, c("CellType", "status", "Time"))))
    for (i in seq_len(nrow(counts))) {
      row <- counts[i, ]
      n <- sum(samples$CellType == row$CellType & samples$status == row$status & samples$Time == row$Time)
      stopifnot(row$N == n, as.logical(row$observed) == (n > 0))
    }
    stopifnot(!any(summary$status == "Died" & summary$Time %in% c("D4", "D8")))
    panels[[panel]] <- list(cells = cells, runs = 184, sample_cell_rows = nrow(samples), mean_points = nrow(summary),
                           no_outlier_filter = TRUE, missing_Died_D4_D8_not_imputed = TRUE)
  }
  a <- get("A/sample_level.tsv")
  b <- get("B/sample_level.tsv")
  same_table(a, b[b$CellType %in% unique(a$CellType), ], c("Run", "CellType"), c("Abundance", "Time", "status"), "A-B shared values")
  early <- b[b$Time %in% c("D0", "D1"), ]
  same_table(get("B/early_model_input.tsv"), early, c("Run", "CellType"), c("Abundance", "Patient_ID", "Time", "TimeDay", "status"), "early model inputs")
  stopifnot(length(unique(early$Run)) == 136, length(unique(early$Patient_ID)) == 135,
            all(table(unique(early[c("Run", "Patient_ID")])$Patient_ID) <= 2),
            sum(table(unique(early[c("Run", "Patient_ID")])$Patient_ID) > 1) == 1)
  source_stats <- get("_shared/B_original_statistics.tsv")
  replay <- list()
  for (window in c("early", "full")) {
    stats <- get(if (window == "early") "B/statistics_displayed.tsv" else "B/statistics_full_window_not_displayed.tsv")
    stopifnot(nrow(stats) == 8, all(stats$Model == "linear_model_lm"), all(stats$Term == "statusSurvived:TimeDay"))
    original <- source_stats[startsWith(source_stats$Analysis, paste0(window, "_")), ]
    same_table(stats, original, "CellType", setdiff(names(stats), "CellType"), paste(window, "source model"))
    values <- if (window == "early") early else b
    recomputed <- stats
    for (i in seq_len(nrow(stats))) {
      dat <- values[values$CellType == stats$CellType[i], ]
      dat$status <- factor(dat$status, c("Died", "Survived"))
      fit <- lm(Abundance ~ status * TimeDay, data = dat)
      coefficient <- summary(fit)$coefficients["statusSurvived:TimeDay", ]
      equal_num(c(stats$Estimate[i], stats$StdError[i], stats$Pvalue[i]),
                unname(coefficient[c(1, 2, 4)]), paste(window, stats$CellType[i], "OLS"))
      stopifnot(stats$N[i] == nrow(dat), stats$N_patient[i] == length(unique(dat$Patient_ID)),
                df.residual(fit) == if (window == "early") 132 else 180)
      recomputed$Estimate[i] <- coefficient[1]
      recomputed$StdError[i] <- coefficient[2]
      recomputed$Pvalue[i] <- coefficient[4]
    }
    recomputed$FDR <- p.adjust(recomputed$Pvalue, "BH")
    equal_num(stats$FDR, recomputed$FDR, paste(window, "eight-test BH FDR"))
    replay[[window]] <- recomputed
  }
  annotations <- get("B/statistics_displayed.tsv")
  stopifnot(identical(sprintf("%.3f", annotations$Pvalue), c("0.002", "0.015", "0.539", "0.321", "0.518", "0.291", "0.022", "0.001")),
            identical(sprintf("%.3f", annotations$FDR), c("0.006", "0.040", "0.539", "0.428", "0.539", "0.428", "0.045", "0.006")))
  list(passed = TRUE, checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), numeric_comparisons = comparisons,
       max_abs_numeric_difference = largest, tolerance = 1e-12, panels = panels,
       frozen_matrix_runs = 211, frozen_state_columns = 53, eligible_runs = 184, early_runs = 136,
       OLS_models_replayed = 16, displayed_P_and_FDR_match_reference_image = TRUE,
       SEM_verified = TRUE, FDR_family = "8 tests per window, not 53", historical_inference_rerun = FALSE,
       recorded_ID_strings_are_not_verified_patient_counts = TRUE, paired_or_mixed_effects_claimed = FALSE,
       original_numeric_values_changed = FALSE, pixel_identity_claimed = FALSE)
}


prepare_analysis <- function(data_dir) {
  list(verification = verify_data(data_dir))
}
