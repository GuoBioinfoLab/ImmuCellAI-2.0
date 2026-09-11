options(stringsAsFactors = FALSE)

fig4_dir <- "<LOCAL_R_ROOT>/Fig4"
cluster_dir <- "<LOCAL_CLUSTER_ROOT>"
workspace_dir <- "<LOCAL_R_ROOT>"

tcga_fraction_file <- file.path(fig4_dir, "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000",
                                "TCGA_ImmunotherapyResponsive_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
tcga_cluster_file <- file.path(cluster_dir, "cluster_TCGA.k=4.consensusClass.csv")
external_expr_file <- file.path(cluster_dir, "RNAseqTpmSymbol.txt")
external_coldata_file <- file.path(cluster_dir, "RNAseq_coldata.Rdata")
reference_file <- file.path(workspace_dir, "reference_53celltypesTPM20260518.txt")
marker_file <- file.path(workspace_dir, "MarkerUsedDeconvolution.txt")

out_dir <- file.path(fig4_dir, "C3_C4_score_external_immunotherapy_validation")
checkpoint_dir <- file.path(out_dir, "checkpoints")
tmp_dir <- file.path(out_dir, "tmp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = tmp_dir, TEMP = tmp_dir, TMPDIR = tmp_dir)

n_cores <- min(8L, parallel::detectCores(logical = TRUE))
chunk_size <- 300L

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ImmuCellAI2.0)
})

log_file <- file.path(out_dir, "validate_C3_C4_score_external_immunotherapy.log")
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

read_gene_tokens <- function(file) {
  txt <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- unlist(strsplit(txt, "[\t,;\"' ]+"))
  x <- trimws(x)
  x <- unique(x[nzchar(x)])
  x <- x[!grepl("^[0-9]+$", x)]
  x <- x[!x %in% c("Gene", "Genes", "gene", "genes", "X", "x", "NA")]
  x
}

read_matrix_gene_first <- function(file, select_samples = NULL) {
  con <- file(file, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1L, warn = FALSE)
  second <- readLines(con, n = 1L, warn = FALSE)
  first_fields <- strsplit(first, "\t", fixed = TRUE)[[1L]]
  second_fields <- strsplit(second, "\t", fixed = TRUE)[[1L]]
  if (length(first_fields) == length(second_fields) - 1L) {
    if (!is.null(select_samples)) {
      keep <- match(select_samples, first_fields)
      keep <- keep[!is.na(keep)]
      dt <- fread(
        file,
        sep = "\t",
        header = FALSE,
        skip = 1L,
        select = c(1L, keep + 1L),
        col.names = c("Gene", first_fields[keep]),
        check.names = FALSE,
        data.table = FALSE
      )
    } else {
      dt <- fread(
        file,
        sep = "\t",
        header = FALSE,
        skip = 1L,
        col.names = c("Gene", first_fields),
        check.names = FALSE,
        data.table = FALSE
      )
    }
  } else {
    dt <- fread(file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
    if (!is.null(select_samples)) {
      cols <- c(colnames(dt)[1], intersect(select_samples, colnames(dt)))
      dt <- dt[, cols, drop = FALSE]
    }
  }
  gene <- as.character(dt[[1L]])
  dt[[1L]] <- NULL
  mat <- as.matrix(dt)
  storage.mode(mat) <- "double"
  rownames(mat) <- gene
  mat[!is.finite(mat)] <- 0
  mat
}

aggregate_duplicate_genes <- function(mat) {
  if (!anyDuplicated(rownames(mat))) return(mat)
  rowsum(mat, group = rownames(mat), reorder = FALSE)
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

cosine_sim <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!length(x) || den == 0) return(NA_real_)
  sum(x * y) / den
}

auc_manual <- function(response_binary, score) {
  ok <- is.finite(score) & !is.na(response_binary)
  y <- response_binary[ok]
  s <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  ranks <- rank(s, ties.method = "average")
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

message("Building TCGA C3/C4 centroids...")
tcga_fraction <- fread(tcga_fraction_file, data.table = FALSE, check.names = FALSE)
rownames(tcga_fraction) <- tcga_fraction[[1L]]
tcga_fraction[[1L]] <- NULL
tcga_fraction <- as.matrix(tcga_fraction)
storage.mode(tcga_fraction) <- "double"

tcga_cluster <- fread(tcga_cluster_file, header = FALSE, data.table = FALSE)
colnames(tcga_cluster) <- c("Sample", "Cluster")
tcga_cluster$Cluster <- paste0("Cluster", as.integer(tcga_cluster$Cluster))
tcga_cluster <- tcga_cluster[tcga_cluster$Sample %in% rownames(tcga_fraction), , drop = FALSE]

tcga_score_mat <- tcga_fraction[tcga_cluster$Sample, , drop = FALSE]
tcga_cluster$Cluster <- factor(tcga_cluster$Cluster, levels = paste0("Cluster", 1:4))
centroids <- sapply(levels(tcga_cluster$Cluster), function(cl) {
  colMeans(tcga_score_mat[tcga_cluster$Cluster == cl, , drop = FALSE], na.rm = TRUE)
})
centroids <- t(centroids)
write.table(centroids, file.path(out_dir, "TCGA_cluster_centroids_53celltypes.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

c3_centroid <- centroids["Cluster3", ]
c4_centroid <- centroids["Cluster4", ]
c4_minus_c3_weight <- c4_centroid - c3_centroid
signature_table <- data.frame(
  CellType = names(c4_minus_c3_weight),
  Cluster3_centroid = as.numeric(c3_centroid),
  Cluster4_centroid = as.numeric(c4_centroid),
  C4_minus_C3 = as.numeric(c4_minus_c3_weight),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(C4_minus_C3))
write.table(signature_table, file.path(out_dir, "TCGA_C4_minus_C3_celltype_signature_weights.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Reading external immunotherapy clinical data...")
load(external_coldata_file)
if (!exists("RNAseq_coldata")) stop("Object RNAseq_coldata not found in: ", external_coldata_file)
coldata <- RNAseq_coldata
coldata$Run <- as.character(coldata$Run)
coldata <- coldata[!is.na(coldata$Run) & nzchar(coldata$Run), , drop = FALSE]

expr_header <- scan(external_expr_file, what = "", nlines = 1L, sep = "\t", quiet = TRUE)
external_samples <- intersect(coldata$Run, expr_header)
if (length(external_samples) < 10L) stop("Too few samples overlap external expression and coldata.")
coldata <- coldata[match(external_samples, coldata$Run), , drop = FALSE]
write.table(coldata, file.path(out_dir, "external_immunotherapy_coldata_matched.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Reading reference and external expression...")
reference <- read_matrix_gene_first(reference_file)
bulk <- read_matrix_gene_first(external_expr_file, select_samples = external_samples)
reference <- aggregate_duplicate_genes(reference)
bulk <- aggregate_duplicate_genes(bulk)

marker_genes <- read_gene_tokens(marker_file)
common <- Reduce(intersect, list(rownames(reference), rownames(bulk), marker_genes))
if (length(common) < 100L) stop("Too few marker genes overlap external bulk and reference: ", length(common))
bulk <- bulk[common, external_samples, drop = FALSE]
reference <- reference[common, , drop = FALSE]
keep <- rowSums(reference) > 0 & rowSums(bulk) > 0
bulk <- bulk[keep, , drop = FALSE]
reference <- reference[keep, , drop = FALSE]

hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = "tcell")
hierarchy <- ImmuCellAI2.0:::validate_hierarchy(hierarchy, reference)
reference <- reference[, hierarchy$state_name, drop = FALSE]

message("Running/reusing ImmuCellAI2 on external immunotherapy samples...")
fraction_file <- file.path(out_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_sample_by_celltype.txt")
if (file.exists(fraction_file)) {
  ext_fraction <- fread(fraction_file, data.table = FALSE, check.names = FALSE)
  rownames(ext_fraction) <- ext_fraction[[1L]]
  ext_fraction[[1L]] <- NULL
  ext_fraction <- as.matrix(ext_fraction)
  storage.mode(ext_fraction) <- "double"
} else {
  samples <- colnames(bulk)
  chunks <- split(seq_along(samples), ceiling(seq_along(samples) / chunk_size))
  state_parts <- vector("list", length(chunks))
  status <- data.frame(
    Chunk = seq_along(chunks),
    Start = vapply(chunks, min, integer(1)),
    End = vapply(chunks, max, integer(1)),
    Samples = vapply(chunks, length, integer(1)),
    Status = "PENDING",
    Seconds = NA_real_
  )
  for (i in seq_along(chunks)) {
    idx <- chunks[[i]]
    checkpoint <- file.path(checkpoint_dir, sprintf("external_chunk_%03d_state_fraction.rds", i))
    if (file.exists(checkpoint)) {
      log_message("Reusing checkpoint chunk ", i, "/", length(chunks))
      state_parts[[i]] <- readRDS(checkpoint)
      status$Status[i] <- "REUSED"
      next
    }
    log_message("Running chunk ", i, "/", length(chunks), " (", length(idx), " samples)...")
    started <- proc.time()[["elapsed"]]
    fit <- deconvolve_bulk_matrix(
      bulk.mat = bulk[, idx, drop = FALSE],
      reference = reference,
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
      seed = 123 + i,
      verbose = FALSE
    )
    part <- as.matrix(fit$state.fraction)
    colnames(part) <- hierarchy$state_name
    state_parts[[i]] <- part
    saveRDS(part, checkpoint, compress = FALSE)
    status$Status[i] <- "OK"
    status$Seconds[i] <- proc.time()[["elapsed"]] - started
    write.table(status, file.path(out_dir, "external_immunotherapy_ImmuCellAI2_chunk_status.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    rm(fit, part)
    gc()
  }
  ext_fraction <- do.call(rbind, state_parts)
  ext_fraction <- ext_fraction[external_samples, hierarchy$state_name, drop = FALSE]
  write.table(ext_fraction, fraction_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(t(ext_fraction), file.path(out_dir, "external_immunotherapy_ImmuCellAI2_marker5000_state_fraction_celltype_by_sample.txt"),
              sep = "\t", quote = FALSE, col.names = NA)
}

message("Scoring external samples...")
common_states <- intersect(colnames(ext_fraction), names(c3_centroid))
score_df <- data.frame(
  Run = rownames(ext_fraction),
  C3_like_spearman = apply(ext_fraction[, common_states, drop = FALSE], 1, safe_cor, y = c3_centroid[common_states], method = "spearman"),
  C4_like_spearman = apply(ext_fraction[, common_states, drop = FALSE], 1, safe_cor, y = c4_centroid[common_states], method = "spearman"),
  C3_like_pearson = apply(ext_fraction[, common_states, drop = FALSE], 1, safe_cor, y = c3_centroid[common_states], method = "pearson"),
  C4_like_pearson = apply(ext_fraction[, common_states, drop = FALSE], 1, safe_cor, y = c4_centroid[common_states], method = "pearson"),
  C3_like_cosine = apply(ext_fraction[, common_states, drop = FALSE], 1, cosine_sim, y = c3_centroid[common_states]),
  C4_like_cosine = apply(ext_fraction[, common_states, drop = FALSE], 1, cosine_sim, y = c4_centroid[common_states]),
  Weighted_C4_minus_C3 = as.numeric(ext_fraction[, common_states, drop = FALSE] %*% c4_minus_c3_weight[common_states]),
  stringsAsFactors = FALSE
)
score_df$C4_minus_C3_spearman <- score_df$C4_like_spearman - score_df$C3_like_spearman
score_df$C4_minus_C3_pearson <- score_df$C4_like_pearson - score_df$C3_like_pearson
score_df$C4_minus_C3_cosine <- score_df$C4_like_cosine - score_df$C3_like_cosine

validation <- coldata %>%
  left_join(score_df, by = "Run") %>%
  mutate(
    ResponseBinary = case_when(
      Response %in% c("R", "Responder", "CR", "PR", "CR/PR") ~ 1L,
      Response %in% c("NR", "Non-responder", "SD", "PD") ~ 0L,
      TRUE ~ NA_integer_
    ),
    ResponseGroup = factor(ifelse(ResponseBinary == 1, "Responder",
                                  ifelse(ResponseBinary == 0, "Non-responder", NA_character_)),
                           levels = c("Non-responder", "Responder"))
  )
write.table(validation, file.path(out_dir, "external_immunotherapy_C3_C4_scores_with_clinical.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

score_vars <- c("C4_minus_C3_spearman", "C4_minus_C3_pearson", "C4_minus_C3_cosine",
                "Weighted_C4_minus_C3", "C4_like_spearman", "C3_like_spearman")

test_one <- function(dat, score_var, group_vars = character()) {
  dat <- dat[!is.na(dat$ResponseBinary) & is.finite(dat[[score_var]]), , drop = FALSE]
  if (nrow(dat) < 6L || length(unique(dat$ResponseBinary)) < 2L) {
    return(data.frame(Score = score_var, N = nrow(dat), N_R = sum(dat$ResponseBinary == 1),
                      N_NR = sum(dat$ResponseBinary == 0), AUC = NA_real_, WilcoxP = NA_real_,
                      stringsAsFactors = FALSE))
  }
  wt <- wilcox.test(dat[[score_var]] ~ dat$ResponseBinary, exact = FALSE)
  data.frame(
    Score = score_var,
    N = nrow(dat),
    N_R = sum(dat$ResponseBinary == 1),
    N_NR = sum(dat$ResponseBinary == 0),
    Mean_R = mean(dat[[score_var]][dat$ResponseBinary == 1], na.rm = TRUE),
    Mean_NR = mean(dat[[score_var]][dat$ResponseBinary == 0], na.rm = TRUE),
    Median_R = median(dat[[score_var]][dat$ResponseBinary == 1], na.rm = TRUE),
    Median_NR = median(dat[[score_var]][dat$ResponseBinary == 0], na.rm = TRUE),
    AUC = auc_manual(dat$ResponseBinary, dat[[score_var]]),
    WilcoxP = wt$p.value,
    stringsAsFactors = FALSE
  )
}

overall_tests <- bind_rows(lapply(score_vars, function(v) test_one(validation, v))) %>%
  mutate(FDR = p.adjust(WilcoxP, method = "BH"))
write.table(overall_tests, file.path(out_dir, "external_immunotherapy_C3_C4_score_response_tests_overall.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

stratified_tests <- bind_rows(lapply(c("Cancer", "Anti_target", "Drug", "SRA_study", "disease"), function(group_var) {
  bind_rows(lapply(split(validation, validation[[group_var]]), function(dat) {
    if (nrow(dat) < 10L || length(unique(na.omit(dat$ResponseBinary))) < 2L) return(NULL)
    res <- bind_rows(lapply(score_vars, function(v) test_one(dat, v)))
    res[[group_var]] <- unique(dat[[group_var]])[1]
    res$GroupVariable <- group_var
    res
  }))
}))
if (nrow(stratified_tests) > 0) {
  stratified_tests <- stratified_tests %>%
    mutate(FDR = p.adjust(WilcoxP, method = "BH")) %>%
    select(GroupVariable, everything())
}
write.table(stratified_tests, file.path(out_dir, "external_immunotherapy_C3_C4_score_response_tests_stratified.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Plotting...")
plot_df <- validation %>%
  filter(!is.na(ResponseGroup), is.finite(C4_minus_C3_spearman))
p_box <- ggplot(plot_df, aes(x = ResponseGroup, y = C4_minus_C3_spearman, fill = ResponseGroup)) +
  geom_boxplot(width = 0.56, outlier.size = 0.6, color = "black", linewidth = 0.25) +
  geom_jitter(width = 0.14, size = 0.45, alpha = 0.45) +
  scale_fill_manual(values = c("Non-responder" = "#4E79A7FF", "Responder" = "#E15759FF"), guide = "none") +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45)) +
  labs(x = NULL, y = "C4-like minus C3-like score")
pdf(file.path(out_dir, "external_immunotherapy_C4_minus_C3_score_response_boxplot.pdf"), width = 4.8, height = 4.5, useDingbats = FALSE)
print(p_box)
dev.off()
png(file.path(out_dir, "external_immunotherapy_C4_minus_C3_score_response_boxplot.png"), width = 1500, height = 1350, res = 300)
print(p_box)
dev.off()

if (requireNamespace("pROC", quietly = TRUE)) {
  roc_obj <- pROC::roc(plot_df$ResponseBinary, plot_df$C4_minus_C3_spearman, quiet = TRUE, direction = "<")
  roc_df <- data.frame(
    Specificity = roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities
  )
  p_roc <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_line(color = "#E15759FF", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    coord_equal() +
    theme_classic(base_size = 12) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45)) +
    labs(x = "1 - Specificity", y = "Sensitivity",
         title = paste0("AUC = ", sprintf("%.3f", as.numeric(pROC::auc(roc_obj)))))
  pdf(file.path(out_dir, "external_immunotherapy_C4_minus_C3_score_response_ROC.pdf"), width = 4.8, height = 4.5, useDingbats = FALSE)
  print(p_roc)
  dev.off()
  png(file.path(out_dir, "external_immunotherapy_C4_minus_C3_score_response_ROC.png"), width = 1500, height = 1350, res = 300)
  print(p_roc)
  dev.off()
}

run_info <- data.frame(
  ExternalSamplesInExpression = length(expr_header),
  ExternalSamplesMatchedClinical = nrow(coldata),
  ExternalSamplesScored = nrow(validation),
  ResponseAvailable = sum(!is.na(validation$ResponseBinary)),
  Responders = sum(validation$ResponseBinary == 1, na.rm = TRUE),
  NonResponders = sum(validation$ResponseBinary == 0, na.rm = TRUE),
  GenesUsed = nrow(bulk),
  States = ncol(ext_fraction),
  Cores = n_cores,
  stringsAsFactors = FALSE
)
write.table(run_info, file.path(out_dir, "external_immunotherapy_C3_C4_validation_run_info.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

log_message("Done. Outputs written to: ", out_dir)
print(run_info)
print(overall_tests)
