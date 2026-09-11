options(stringsAsFactors = FALSE)

workspace_dir <- "<LOCAL_R_ROOT>"
fig4_dir <- file.path(workspace_dir, "Fig4")
out_dir <- file.path(fig4_dir, "gene_set_performance_compare")
tmp_dir <- file.path(workspace_dir, "Rtmp_immucellai2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = tmp_dir, TEMP = tmp_dir, TMPDIR = tmp_dir)

reference_file <- file.path(workspace_dir, "reference_53celltypesTPM20260518.txt")

gene_set_files <- c(
  Immune_5499_gene_set = file.path(workspace_dir, "MarkerUsedDeconvolution.txt"),
  Fig4_marker_gene_set = file.path(fig4_dir, "MarkerUsedDeconvolution.txt")
)

datasets <- data.frame(
  Dataset = c("GSE146771", "GSE164522", "GSE176078"),
  BulkFile = c(
    "<LOCAL_PROJECT_ROOT>/GSE146771_simulate_variable_background.txt",
    "<LOCAL_PROJECT_ROOT>/GSE164522_simulate_variable_background.txt",
    "<LOCAL_PROJECT_ROOT>/GSE176078_simulate_variable_background.txt"
  ),
  stringsAsFactors = FALSE
)

n_cores <- min(8L, parallel::detectCores(logical = TRUE))

suppressPackageStartupMessages({
  library(data.table)
  library(ImmuCellAI2.0)
})

read_expression_matrix <- function(file) {
  con <- file(file, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1L, warn = FALSE)
  second <- readLines(con, n = 1L, warn = FALSE)
  first_fields <- strsplit(first, "\t", fixed = TRUE)[[1L]]
  second_fields <- strsplit(second, "\t", fixed = TRUE)[[1L]]

  if (length(first_fields) == length(second_fields) - 1L) {
    x <- data.table::fread(
      file,
      sep = "\t",
      header = FALSE,
      skip = 1L,
      col.names = c("Gene", first_fields),
      check.names = FALSE,
      data.table = FALSE
    )
  } else {
    x <- data.table::fread(file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  }
  gene <- as.character(x[[1L]])
  x[[1L]] <- NULL
  mat <- as.matrix(x)
  storage.mode(mat) <- "double"
  rownames(mat) <- gene
  mat[!is.finite(mat)] <- 0
  mat
}

read_gene_tokens <- function(file) {
  if (!file.exists(file)) stop("Gene-set file not found: ", file)
  txt <- readLines(file, warn = FALSE, encoding = "UTF-8")
  tokens <- unlist(strsplit(txt, "[\t,;\"' ]+"))
  tokens <- trimws(tokens)
  tokens <- unique(tokens[nzchar(tokens)])
  tokens <- tokens[!grepl("^[0-9]+$", tokens)]
  tokens <- tokens[!tokens %in% c("Gene", "Genes", "gene", "genes", "X", "x", "NA")]
  tokens
}

aggregate_duplicate_genes <- function(mat) {
  if (!anyDuplicated(rownames(mat))) return(mat)
  rowsum(mat, group = rownames(mat), reorder = FALSE)
}

extract_variable_sample_info <- function(samples) {
  parsed <- regexec("^(.+)_([0-9]+)$", samples)
  pieces <- regmatches(samples, parsed)
  ok <- lengths(pieces) == 3L
  if (!all(ok)) {
    stop("Cannot parse sample names: ", paste(utils::head(samples[!ok], 10), collapse = ", "))
  }
  data.frame(
    Sample = samples,
    CellType = vapply(pieces, `[`, character(1), 2),
    StepIndex = as.integer(vapply(pieces, `[`, character(1), 3)),
    TrueFraction = as.integer(vapply(pieces, `[`, character(1), 3)) * 0.005,
    stringsAsFactors = FALSE
  )
}

state_target_mapper <- function(target) {
  mapper <- list(
    Bcells = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC"),
    Bmemory = c("MBC", "Bex"),
    CD4Tmemory = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm"),
    CD4Tcm = "CD4Tcm",
    CD4Tem = "CD4Tem",
    CD4Temra = "CD4Temra",
    CD4Tfh = "Tfh",
    CD4Tnaive = "CD4Tn",
    CD4Trm = "CD4Trm",
    CD8Tcm = "CD8Tcm",
    CD8Tem = "CD8Tem",
    CD8Temra = "CD8Temra",
    CD8Tnaive = "CD8Tn",
    CD8Trm = "CD8Trm",
    CD8Tmemory = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm"),
    CMonocyte = "cMo",
    GC_B = "BGC",
    ILC3 = "ILC3",
    MAIT = "MAIT",
    Monocyte = c("cMo", "intMo", "ncMo"),
    NMonocyte = "ncMo",
    NK = c("cNK", "NKreg"),
    NKT = "NKT",
    TAM = "TAM",
    Thelper = c("Tfh", "Th1", "Th1/17", "Th17", "Th2"),
    Th1 = "Th1",
    Th17 = "Th17",
    Treg = "Treg",
    cDC1 = "cDC1",
    cDC2 = "cDC2",
    cytotoxicNK = "cNK",
    exhausted_T = "Tex",
    macrophage = c("M0", "M1", "M2"),
    mast_cell = "Mast cell",
    monoDC = "monoDC",
    pDC = "pDC",
    plasma = c("PB", "PC"),
    regulatoryNK = "NKreg",
    gdT = "gdT",
    `γδT` = "gdT"
  )
  out <- lapply(as.character(target), function(x) if (x %in% names(mapper)) mapper[[x]] else x)
  names(out) <- target
  out
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

evaluate_prediction <- function(pred, sample_info) {
  pred <- as.matrix(pred)
  if (!all(sample_info$Sample %in% rownames(pred)) && all(sample_info$Sample %in% colnames(pred))) {
    pred <- t(pred)
  }
  pred <- pred[sample_info$Sample, , drop = FALSE]
  targets <- unique(sample_info$CellType)
  rows <- lapply(targets, function(ct) {
    states <- intersect(state_target_mapper(ct)[[ct]], colnames(pred))
    idx <- sample_info$CellType == ct
    truth <- sample_info$TrueFraction[idx]
    value <- if (length(states) == 0L) rep(NA_real_, sum(idx)) else rowSums(pred[idx, states, drop = FALSE], na.rm = TRUE)
    ok <- is.finite(value) & is.finite(truth)
    fit <- if (sum(ok) >= 3L && stats::sd(value[ok]) > 0 && stats::sd(truth[ok]) > 0) {
      stats::coef(stats::lm(value[ok] ~ truth[ok]))
    } else {
      c(`(Intercept)` = NA_real_, truth = NA_real_)
    }
    data.frame(
      CellType = ct,
      PredictedStates = paste(states, collapse = "+"),
      Pearson = safe_cor(value, truth, "pearson"),
      Spearman = safe_cor(value, truth, "spearman"),
      RMSE = if (sum(ok) > 0L) sqrt(mean((value[ok] - truth[ok])^2)) else NA_real_,
      MAE = if (sum(ok) > 0L) mean(abs(value[ok] - truth[ok])) else NA_real_,
      Slope = unname(fit[2]),
      Intercept = unname(fit[1]),
      N = sum(ok),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

standardise_fraction_matrix <- function(pred, state_names, sample_names = NULL) {
  pred <- as.matrix(pred)
  storage.mode(pred) <- "double"
  if (ncol(pred) == length(state_names)) {
    colnames(pred) <- state_names
  } else if (nrow(pred) == length(state_names)) {
    rownames(pred) <- state_names
  }
  if (!is.null(sample_names) && nrow(pred) == length(sample_names) && is.null(rownames(pred))) {
    rownames(pred) <- sample_names
  }
  pred
}

read_fraction_matrix <- function(file, state_names, sample_names) {
  raw <- data.table::fread(file, sep = "\t", header = FALSE, check.names = FALSE, data.table = FALSE)
  first_col <- as.character(raw[[1L]])
  if (length(first_col) > 0L && !nzchar(first_col[1L])) {
    raw <- raw[-1L, , drop = FALSE]
    first_col <- first_col[-1L]
  }
  rn <- as.character(raw[[1L]])
  raw[[1L]] <- NULL
  pred <- as.matrix(raw)
  storage.mode(pred) <- "double"
  rownames(pred) <- rn
  standardise_fraction_matrix(pred, state_names = state_names, sample_names = sample_names)
}

summarise_metrics <- function(metrics, dataset, gene_set, genes_used, seconds) {
  data.frame(
    Dataset = dataset,
    GeneSet = gene_set,
    GenesUsed = genes_used,
    Seconds = seconds,
    NTargets = nrow(metrics),
    ValidPearsonTargets = sum(!is.na(metrics$Pearson)),
    MeanPearson = mean(metrics$Pearson, na.rm = TRUE),
    MedianPearson = stats::median(metrics$Pearson, na.rm = TRUE),
    MeanSpearman = mean(metrics$Spearman, na.rm = TRUE),
    MedianSpearman = stats::median(metrics$Spearman, na.rm = TRUE),
    MeanRMSE = mean(metrics$RMSE, na.rm = TRUE),
    MedianRMSE = stats::median(metrics$RMSE, na.rm = TRUE),
    MeanMAE = mean(metrics$MAE, na.rm = TRUE),
    MedianMAE = stats::median(metrics$MAE, na.rm = TRUE),
    MeanSlope = mean(metrics$Slope, na.rm = TRUE),
    MedianSlope = stats::median(metrics$Slope, na.rm = TRUE),
    MeanIntercept = mean(metrics$Intercept, na.rm = TRUE),
    MedianIntercept = stats::median(metrics$Intercept, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

plot_summary <- function(summary_table, file) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  pdat <- summary_table
  pdat$GeneSet <- factor(pdat$GeneSet, levels = names(gene_set_files))
  p <- ggplot2::ggplot(pdat, ggplot2::aes(x = Dataset, y = MeanPearson, fill = GeneSet)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", MeanPearson)),
      position = ggplot2::position_dodge(width = 0.75),
      vjust = -0.25,
      size = 3.2
    ) +
    ggplot2::coord_cartesian(ylim = c(max(0, min(pdat$MeanPearson, na.rm = TRUE) - 0.04), 1.02)) +
    ggplot2::scale_fill_manual(values = c("#F66463", "#379DA5")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.45),
      legend.position = "top"
    ) +
    ggplot2::labs(x = NULL, y = "Mean Pearson", fill = NULL)
  ggplot2::ggsave(file, p, width = 7.2, height = 4.2)
  invisible(TRUE)
}

message("Reading reference: ", reference_file)
reference_all <- read_expression_matrix(reference_file)
reference_all <- aggregate_duplicate_genes(reference_all)

gene_sets <- lapply(gene_set_files, read_gene_tokens)
gene_set_info <- data.frame(
  GeneSet = names(gene_sets),
  RawTokens = vapply(gene_sets, length, integer(1)),
  ReferenceOverlap = vapply(gene_sets, function(x) length(intersect(x, rownames(reference_all))), integer(1)),
  stringsAsFactors = FALSE
)
write.table(gene_set_info, file.path(out_dir, "gene_set_overlap_info.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(gene_set_info)

all_metrics <- list()
all_summary <- list()

for (d in seq_len(nrow(datasets))) {
  dataset <- datasets$Dataset[d]
  bulk_file <- datasets$BulkFile[d]
  if (!file.exists(bulk_file)) stop("Bulk file not found: ", bulk_file)

  message("\n=== ", dataset, " ===")
  bulk_all <- read_expression_matrix(bulk_file)
  bulk_all <- aggregate_duplicate_genes(bulk_all)
  sample_info <- extract_variable_sample_info(colnames(bulk_all))

  common <- intersect(rownames(bulk_all), rownames(reference_all))
  bulk_common <- bulk_all[common, , drop = FALSE]
  ref_common <- reference_all[common, , drop = FALSE]
  keep <- rowSums(bulk_common, na.rm = TRUE) > 0 & rowSums(ref_common, na.rm = TRUE) > 0
  bulk_common <- bulk_common[keep, , drop = FALSE]
  ref_common <- ref_common[keep, , drop = FALSE]
  rm(bulk_all)
  gc()

  hierarchy <- create_default_53_hierarchy(colnames(ref_common), hierarchy.mode = "tcell")
  hierarchy <- ImmuCellAI2.0:::validate_hierarchy(hierarchy, ref_common)
  ref_common <- ref_common[, hierarchy$state_name, drop = FALSE]

  for (gene_set_name in names(gene_sets)) {
    candidate_genes <- intersect(gene_sets[[gene_set_name]], rownames(ref_common))
    if (length(candidate_genes) < 100L) {
      warning("Too few genes for ", dataset, " / ", gene_set_name, ": ", length(candidate_genes))
      next
    }

    fraction_file <- file.path(out_dir, paste0(dataset, "_", gene_set_name, "_state_fraction.txt"))
    metrics_file <- file.path(out_dir, paste0(dataset, "_", gene_set_name, "_per_celltype_metrics.txt"))

    message("Running/reusing ", gene_set_name, " with ", length(candidate_genes), " genes")
    started <- proc.time()[["elapsed"]]
    if (file.exists(fraction_file)) {
      pred <- read_fraction_matrix(
        fraction_file,
        state_names = hierarchy$state_name,
        sample_names = colnames(bulk_common)
      )
      write.table(pred, fraction_file, sep = "\t", quote = FALSE, col.names = NA)
      seconds <- NA_real_
    } else {
      fit <- deconvolve_bulk_matrix(
        bulk.mat = bulk_common[candidate_genes, , drop = FALSE],
        reference = ref_common[candidate_genes, , drop = FALSE],
        hierarchy = hierarchy,
        hierarchy.mode = "tcell",
        inference.method = "vb",
        add.unknown = FALSE,
        pseudo.depth = 1e5,
        n.iter = 50,
        vb.tol = 1e-6,
        alpha.major = 10,
        alpha.sub = 5,
        alpha.state = 1,
        n.cores = n_cores,
        seed = 123,
        verbose = FALSE
      )
      pred <- standardise_fraction_matrix(
        fit$state.fraction,
        state_names = hierarchy$state_name,
        sample_names = colnames(bulk_common)
      )
      seconds <- proc.time()[["elapsed"]] - started
      write.table(pred, fraction_file, sep = "\t", quote = FALSE, col.names = NA)
      rm(fit)
      gc()
    }

    metrics <- evaluate_prediction(pred, sample_info)
    metrics$Dataset <- dataset
    metrics$GeneSet <- gene_set_name
    metrics$GenesUsed <- length(candidate_genes)
    metrics <- metrics[, c("Dataset", "GeneSet", "GenesUsed", setdiff(colnames(metrics), c("Dataset", "GeneSet", "GenesUsed")))]
    write.table(metrics, metrics_file, sep = "\t", quote = FALSE, row.names = FALSE)

    key <- paste(dataset, gene_set_name, sep = "__")
    all_metrics[[key]] <- metrics
    all_summary[[key]] <- summarise_metrics(metrics, dataset, gene_set_name, length(candidate_genes), seconds)
  }

  rm(bulk_common, ref_common)
  gc()
}

combined_metrics <- do.call(rbind, all_metrics)
combined_summary <- do.call(rbind, all_summary)
rownames(combined_metrics) <- NULL
rownames(combined_summary) <- NULL
combined_summary <- combined_summary[order(combined_summary$Dataset, combined_summary$GeneSet), ]

write.table(combined_metrics, file.path(out_dir, "all_datasets_gene_set_per_celltype_metrics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(combined_summary, file.path(out_dir, "all_datasets_gene_set_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

wide <- reshape(
  combined_metrics[, c("Dataset", "CellType", "GeneSet", "Pearson")],
  idvar = c("Dataset", "CellType"),
  timevar = "GeneSet",
  direction = "wide"
)
colnames(wide) <- sub("^Pearson\\.", "", colnames(wide))
write.table(wide, file.path(out_dir, "all_datasets_gene_set_pearson_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_summary(combined_summary, file.path(out_dir, "gene_set_mean_pearson_compare.pdf"))

message("\nDone. Results written to: ", out_dir)
print(combined_summary)
