analysis$E <- local({
  get <- function(path) read_tsv(file.path(data_dir, path))
  metadata <- inputs$metadata
  check_values <- inputs$check_values
  grouping <- get("_shared/original_metadata/grouping_sheet1.tsv")
  stopifnot(!anyDuplicated(grouping$Run), setequal(grouping$Run, metadata$PRJNA683803$Run))
  mapping <- get("E/time_mapping.tsv")
  stopifnot(!anyDuplicated(mapping$grouping))
  meta <- metadata$PRJNA683803
  meta$grouping <- grouping$grouping[match(meta$Run, grouping$Run)]
  meta$Time <- mapping$Time[match(meta$grouping, mapping$grouping)]
  eligible <- meta[meta$outcome %in% c("Died", "Survived") & !is.na(meta$Time), ]
  samples <- get("E/sample_level.tsv")
  cells <- get("E/cell_order.tsv")$CellType
  stopifnot(nrow(samples) == nrow(eligible) * length(cells), setequal(samples$Run, eligible$Run),
            !anyDuplicated(key(samples, c("Run", "CellType"))), identical(unique(samples$CellType), cells))
  check_values(samples, "PRJNA683803", "outcome")
  stopifnot(identical(samples$Time, eligible$Time[match(samples$Run, eligible$Run)]))
  means <- get("E/plot_data.tsv")
  expected <- aggregate(Abundance ~ CellType + status + Time, samples, mean)
  names(expected)[4] <- "Mean"
  check_table(means, expected, c("CellType", "status", "Time"), "Mean", "E means")
  check_table(means, get("E/original_group_means.tsv"), c("CellType", "status", "Time"), "Mean", "E archived means")
  stopifnot(nrow(means) == 24, !any(means$status == "Died" & means$Time %in% c("D4", "D8")))
  counts <- get("E/sample_counts.tsv")
  for (i in seq_len(nrow(counts))) {
    row <- counts[i, ]
    n <- sum(samples$CellType == row$CellType & samples$status == row$status & samples$Time == row$Time)
    stopifnot(row$n == n, (tolower(as.character(row$observed)) == "true") == (n > 0))
  }
  report <- list(sample_cell_rows = nrow(samples), input_runs = nrow(eligible), mean_points = nrow(means),
                         no_outlier_filter = TRUE, died_D4_D8_absent_not_zero = TRUE)

  list(plot_data = means, cells = cells, groups = get("E/group_order.tsv")$status, report = report)
})
