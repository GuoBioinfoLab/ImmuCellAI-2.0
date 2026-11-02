read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             fileEncoding = "UTF-8-BOM", quote = "", comment.char = "")
}

key <- function(dat, cols) do.call(paste, c(dat[cols], sep = "|"))

comparisons <- 0L
largest_difference <- 0
check_num <- function(actual, expected, label, tolerance = 1e-12) {
  stopifnot(length(actual) == length(expected), identical(is.na(actual), is.na(expected)))
  difference <- if (length(actual)) max(abs(actual - expected), na.rm = TRUE) else 0
  if (!is.finite(difference)) difference <- 0
  if (difference > tolerance) stop(label, ": numerical difference ", difference)
  largest_difference <<- max(largest_difference, difference)
  comparisons <<- comparisons + 1L
}
check_table <- function(actual, expected, keys, numeric_cols, label) {
  ka <- key(actual, keys)
  ke <- key(expected, keys)
  stopifnot(!anyDuplicated(ka), !anyDuplicated(ke), setequal(ka, ke))
  expected <- expected[match(ka, ke), , drop = FALSE]
  for (name in numeric_cols) check_num(actual[[name]], expected[[name]], paste(label, name))
}

verify_inputs <- function(data_dir) {
  get <- function(path) read_tsv(file.path(data_dir, path))
  states <- metadata <- list()
  for (cohort in c("TB", "PRJNA683803")) {
    dat <- get(paste0("_shared/", cohort, "/state_fraction.tsv"))
    names(dat)[1] <- "Run"
    meta <- get(paste0("_shared/", cohort, "/metadata.tsv"))
    stopifnot(ncol(dat) == 54, !anyDuplicated(dat$Run), !anyDuplicated(meta$Run),
              setequal(dat$Run, meta$Run), nrow(dat) == if (cohort == "TB") 92 else 211)
    values <- as.matrix(dat[-1])
    stopifnot(all(is.finite(values)), all(values >= 0), all(values <= 1))
    check_num(rowSums(values), rep(1, nrow(dat)), paste(cohort, "fraction sums"))
    states[[cohort]] <- dat
    metadata[[cohort]] <- meta
  }
  check_values <- function(dat, cohort, group_column = NULL) {
    state <- states[[cohort]]
    meta <- metadata[[cohort]]
    stopifnot(all(dat$Run %in% state$Run), all(dat$CellType %in% names(state)))
    expected <- vapply(seq_len(nrow(dat)), function(i) state[[dat$CellType[i]]][match(dat$Run[i], state$Run)], numeric(1))
    check_num(dat$Abundance, expected, paste(cohort, "sample-level values"))
    if (!is.null(group_column)) {
      stopifnot(identical(dat$status, meta[[group_column]][match(dat$Run, meta$Run)]))
    }
  }
  original_tb <- read.csv(file.path(data_dir, "_shared/original_metadata/PRJNA_combined.csv"), check.names = FALSE)
  original_aids <- read.csv(file.path(data_dir, "_shared/original_metadata/PRJNA683803_combined.csv"), check.names = FALSE)
  for (cohort in names(metadata)) {
    original <- if (cohort == "TB") original_tb else original_aids
    original <- original[!duplicated(original$Run), , drop = FALSE]
    m <- metadata[[cohort]]
    stopifnot(all(m$Run %in% original$Run))
    original <- original[match(m$Run, original$Run), , drop = FALSE]
    for (name in names(m)) stopifnot(identical(as.character(m[[name]]), as.character(original[[name]])))
  }
  combined <- get("_shared/PRJNA683803/state_fraction_with_metadata.tsv")
  state <- states$PRJNA683803
  stopifnot(!anyDuplicated(combined$Run), setequal(combined$Run, state$Run))
  for (name in names(state)[-1]) check_num(combined[[name]], state[[name]][match(combined$Run, state$Run)], name)

  is_outlier <- function(x) {
    result <- rep(FALSE, length(x))
    ok <- is.finite(x)
    if (sum(ok) < 4 || length(unique(x[ok])) < 2) return(result)
    quartiles <- quantile(x[ok], c(0.25, 0.75), names = FALSE, type = 7)
    spread <- diff(quartiles)
    if (!is.finite(spread) || spread <= 0) return(result)
    result[ok] <- x[ok] < quartiles[1] - 1.5 * spread | x[ok] > quartiles[2] + 1.5 * spread
    result
  }
  symbol <- function(p) ifelse(is.na(p), "ns", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))
  all_stats <- filtered <- raw_tables <- list()
  recomputed_stats <- list()
  for (kind in c("outcome", "treatment")) {
    raw <- get(paste0("_shared/PRJNA683803/", kind, "_all53_raw.tsv"))
    kept <- get(paste0("_shared/PRJNA683803/", kind, "_all53_filtered.tsv"))
    stopifnot(nrow(raw) == 211 * 53, !anyDuplicated(key(raw, c("Run", "CellType"))))
    check_values(raw, "PRJNA683803", kind)
    flags <- rep(FALSE, nrow(raw))
    for (rows in split(seq_len(nrow(raw)), key(raw, c("CellType", "status")))) flags[rows] <- is_outlier(raw$Abundance[rows])
    check_table(kept, raw[!flags, ], c("Run", "CellType", "status"), "Abundance", paste(kind, "filter"))
    stopifnot(all(!kept$IsOutlier))
    stats <- get(paste0("_shared/PRJNA683803/", kind, "_all53_stats.tsv"))
    stopifnot(nrow(stats) == if (kind == "outcome") 106 else 53,
              !anyDuplicated(key(stats, c("CellType", "group1", "group2"))))
    recomputed <- stats
    for (i in seq_len(nrow(stats))) {
      row <- stats[i, ]
      x <- kept$Abundance[kept$CellType == row$CellType & kept$status == row$group1]
      y <- kept$Abundance[kept$CellType == row$CellType & kept$status == row$group2]
      p <- suppressWarnings(wilcox.test(x, y)$p.value)
      expected <- c(length(x), length(y), mean(x), mean(y), median(x), median(y),
                    log2((mean(y) + 1e-12) / (mean(x) + 1e-12)), p)
      fields <- c("n1", "n2", "mean1", "mean2", "median1", "median2", "log2FC_group2_vs_group1", "p")
      check_num(as.numeric(row[fields]), expected, paste(kind, row$CellType, row$group1))
      stopifnot(row$p.signif == symbol(p))
      recomputed$p[i] <- p
    }
    for (rows in split(seq_len(nrow(stats)), key(stats, c("group1", "group2")))) {
      stopifnot(length(rows) == 53)
      check_num(stats$FDR[rows], p.adjust(recomputed$p[rows], "BH"), paste(kind, "53-cell BH FDR"))
    }
    recomputed_stats[[kind]] <- recomputed
    all_stats[[kind]] <- stats
    filtered[[kind]] <- kept
    raw_tables[[kind]] <- raw
  }
  list(states = states, metadata = metadata, filtered = filtered,
       raw_tables = raw_tables, all_stats = all_stats,
       recomputed_stats = recomputed_stats, check_values = check_values)
}

analyze_box_panel <- function(panel, inputs, data_dir) {
  get <- function(path) read_tsv(file.path(data_dir, path))
  check_values <- inputs$check_values
  raw_tables <- inputs$raw_tables
  filtered <- inputs$filtered
  all_stats <- inputs$all_stats
  symbol <- function(p) ifelse(is.na(p), "ns", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns"))))
  index <- get("panel_index.tsv")
  stopifnot(identical(index$panel, LETTERS[1:5]))
  cells <- get(paste0(panel, "/cell_order.tsv"))$CellType
  groups <- get(paste0(panel, "/group_order.tsv"))$status
  plotted <- get(paste0(panel, "/plot_data.tsv"))
  raw <- get(paste0(panel, "/all_values_with_display_flags.tsv"))
  stats <- get(paste0(panel, "/statistics.tsv"))
  stopifnot(identical(unique(plotted$CellType), cells), identical(stats$CellType, cells))
  check_values(raw, if (panel == "A") "TB" else "PRJNA683803",
               if (panel == "A") "tb_status" else if (panel == "D") "treatment" else "outcome")
  if (panel == "A") {
    stopifnot(nrow(raw) == 368, all(tolower(as.character(raw$IncludedInPlot)) == "true"), all(!plotted$IsPlotOutlier))
    check_table(plotted, raw, c("Run", "CellType", "status"), "Abundance", "A all values")
    for (i in seq_len(nrow(stats))) {
      x <- plotted$Abundance[plotted$CellType == cells[i] & plotted$status == groups[1]]
      y <- plotted$Abundance[plotted$CellType == cells[i] & plotted$status == groups[2]]
      check_num(as.numeric(stats[i, c("N_Active", "N_Latent", "Mean_Active", "Mean_Latent", "P_Wilcox")]),
                c(length(x), length(y), mean(x), mean(y), wilcox.test(x, y)$p.value), paste("A stats", cells[i]))
    }
    expected_stars <- rep("***", 4)
    stopifnot(identical(symbol(stats$P_Wilcox), expected_stars))
  } else {
    kind <- if (panel == "D") "treatment" else "outcome"
    select_rows <- function(dat) dat[dat$CellType %in% cells & dat$status %in% groups, , drop = FALSE]
    check_table(raw, select_rows(raw_tables[[kind]]), c("Run", "CellType", "status"), "Abundance", paste(panel, "raw"))
    check_table(plotted, select_rows(filtered[[kind]]), c("Run", "CellType", "status"), "Abundance", paste(panel, "filtered"))
    expected_stats <- all_stats[[kind]]
    expected_stats <- expected_stats[expected_stats$CellType %in% cells & expected_stats$group1 == groups[1] & expected_stats$group2 == groups[2], ]
    numeric_cols <- names(stats)[vapply(stats, is.numeric, logical(1))]
    check_table(stats, expected_stats, c("CellType", "group1", "group2"), numeric_cols, paste(panel, "stats subset"))
    expected_stars <- if (panel == "C") c("ns", "ns", "*", "ns") else rep("***", 4)
    stopifnot(identical(stats$p.signif, expected_stars))
  }
  included <- tolower(as.character(raw$IncludedInPlot)) == "true"
  stopifnot(identical(included, key(raw, c("Run", "CellType")) %in% key(plotted, c("Run", "CellType"))))
  counts <- get(paste0(panel, "/sample_counts.tsv"))
  for (i in seq_len(nrow(counts))) {
    count <- counts[i, ]
    n_raw <- sum(raw$CellType == count$CellType & raw$status == count$status)
    n_plot <- sum(plotted$CellType == count$CellType & plotted$status == count$status)
    stopifnot(count$n_raw == n_raw, count$n_plotted == n_plot, count$n_excluded == n_raw - n_plot)
  }
  stopifnot(nrow(plotted) == index$plot_rows[index$panel == panel])
  report <- list(raw = nrow(raw), plotted = nrow(plotted), excluded = nrow(raw) - nrow(plotted),
                               input_runs = length(unique(raw$Run)), stars = expected_stars)

  list(plot_data = plotted, statistics = stats, cells = cells, groups = groups, report = report)
}
