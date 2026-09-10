# Imported original Figure 3 analysis; see ../../README.md for provenance.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

options(stringsAsFactors = FALSE)

workspace_dir <- fig3_work_dir()
runtime_compare_dir <- file.path(workspace_dir, "Fig2", "runtime_compare")
local_lib <- file.path(workspace_dir, "Rlib_test")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

stable_tmp <- file.path(workspace_dir, "Rtmp_immucellai2")
dir.create(stable_tmp, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMP = stable_tmp, TEMP = stable_tmp, TMPDIR = stable_tmp)

reference_file <- fig3_reference()
immune_gene_file <- fig3_markers()

env_samples <- Sys.getenv("IMMUCELLAI2_RUNTIME_SAMPLES", "")
env_cores <- Sys.getenv("IMMUCELLAI2_RUNTIME_CORES", "")
n_benchmark_samples <- if (nzchar(env_samples)) as.integer(env_samples) else 12
n_cores <- if (nzchar(env_cores)) {
  as.integer(env_cores)
} else {
  min(8, parallel::detectCores(logical = TRUE))
}
n_benchmark_samples <- max(1, n_benchmark_samples)
n_cores <- max(1, min(n_cores, parallel::detectCores(logical = TRUE)))
scenario_label <- Sys.getenv(
  "IMMUCELLAI2_RUNTIME_LABEL",
  paste0("n", n_benchmark_samples, "_cores", n_cores)
)
out_dir <- file.path(runtime_compare_dir, paste0("runtime_benchmark_7tools_", scenario_label))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seed <- 20260604
chain_length <- 600
burn_in <- 500
cibersort_signature_genes <- 1000

need_pkgs <- c("ImmuCellAI2.0", "BayesPrism", "CITMIC", "ImmuCellAI",
               "MuSiC", "e1071", "limma", "ggplot2")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

library(ImmuCellAI2.0)
library(BayesPrism)
library(CITMIC)
library(ImmuCellAI)
library(MuSiC)
library(ggplot2)

patch_immucellai_gsva <- function() {
  ns <- asNamespace("ImmuCellAI")
  if (!exists("Sample_abundance_calculation", envir = ns, inherits = FALSE)) {
    return(invisible(FALSE))
  }
  patched_sample_abundance <- function(sample, paper_marker, marker_exp, data_type, customer) {
    sam <- apply(sample, 2, as.numeric)
    row.names(sam) <- row.names(sample)
    tt <- intersect(row.names(sam), as.vector(unlist(paper_marker)))
    genes <- intersect(tt, row.names(marker_exp))
    sam_exp <- as.matrix(sam[genes, , drop = FALSE])
    colnames(sam_exp) <- colnames(sample)
    marker_exp <- marker_exp[genes, , drop = FALSE]
    marker_tag_mat <- c()
    for (cell in names(paper_marker)) {
      tag <- marker_tag(genes, as.vector(unlist(paper_marker[cell])))
      marker_tag_mat <- cbind(marker_tag_mat, tag)
    }
    row.names(marker_tag_mat) <- row.names(marker_exp)
    colnames(marker_tag_mat) <- names(paper_marker)
    exp_new <- apply(sam_exp, 2, function(x) sample_ratio(x, marker_exp, marker_tag_mat, data_type))
    if (requireNamespace("GSVA", quietly = TRUE) &&
        exists("ssgseaParam", envir = asNamespace("GSVA"))) {
      param <- GSVA::ssgseaParam(exprData = as.matrix(exp_new),
                                 geneSets = paper_marker,
                                 normalize = TRUE,
                                 verbose = FALSE)
      result <- GSVA::gsva(param)
    } else {
      result <- GSVA::gsva(as.matrix(exp_new), paper_marker, method = "ssgsea",
                           ssgsea.norm = TRUE, parallel.sz = min(20, n_cores))
    }
    if (ncol(result) < 3) {
      result[which(result < 0)] <- 0
    } else {
      result <- result - apply(result, 1, min)
    }
    compensation_matrix_num <- apply(compensation_matrix, 2, as.numeric)
    if (customer == 0) {
      row.names(compensation_matrix_num) <- row.names(compensation_matrix)
      result_norm <- compensation(result, compensation_matrix_num)
      if (ncol(result_norm) == 1) {
        InfiltrationScore <- sum(result_norm)
      } else if (nrow(result_norm) == 24) {
        InfiltrationScore <- apply(
          result_norm[c("Bcell", "CD4_T", "CD8_T", "DC", "Macrophage", "Monocyte", "Neutrophil", "NK"), ],
          2, sum
        )
      } else {
        InfiltrationScore <- apply(result_norm, 2, sum)
      }
      InfiltrationScore <- (InfiltrationScore / max(InfiltrationScore)) * 0.9
      result_mat <- rbind(result_norm, InfiltrationScore)
    } else {
      result_mat <- result
    }
    if (ncol(result_mat) > 1) {
      result_mat <- apply(result_mat, 1, function(x) round(x, 3))
    } else {
      result_mat <- t(round(result_mat, 3))
    }
    result_mat <- as.matrix(result_mat)
    if (ncol(result_mat) == 1) {
      colnames(result_mat) <- colnames(sample)
      t(result_mat)
    } else {
      result_mat
    }
  }
  environment(patched_sample_abundance) <- ns
  was_locked <- bindingIsLocked("Sample_abundance_calculation", ns)
  if (was_locked) unlockBinding("Sample_abundance_calculation", ns)
  assign("Sample_abundance_calculation", patched_sample_abundance, envir = ns)
  if (was_locked) lockBinding("Sample_abundance_calculation", ns)
  invisible(TRUE)
}

read_expression_matrix <- function(file) {
  x <- read.delim(file, row.names = 1, check.names = FALSE, sep = "\t")
  mat <- as.matrix(x)
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat
}

ensure_r_tempdir <- function() {
  td <- tempdir()
  if (!dir.exists(td)) {
    dir.create(td, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(td)
}

read_gene_list <- function(file) {
  if (!file.exists(file)) return(character(0))
  x <- readLines(file, warn = FALSE, encoding = "UTF-8")
  x <- unlist(strsplit(x, "[\t,; ]+"))
  x <- trimws(x)
  unique(x[nzchar(x)])
}

normalize_cols <- function(mat, eps = 1e-12) {
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  cs <- colSums(mat)
  cs[!is.finite(cs) | cs <= 0] <- eps
  sweep(mat, 2, cs, "/")
}

normalize_rows <- function(mat, eps = 1e-12) {
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  rs <- rowSums(mat)
  rs[!is.finite(rs) | rs <= 0] <- eps
  sweep(mat, 1, rs, "/")
}

normalize_nonnegative <- function(x) {
  x[!is.finite(x)] <- 0
  x <- pmax(x, 0)
  s <- sum(x)
  if (s <= 0) rep(1 / length(x), length(x)) else x / s
}

weighted_lm_nonnegative <- function(A, b, w = NULL) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  b[!is.finite(b)] <- 0
  if (is.null(w)) {
    Aw <- A
    bw <- b
  } else {
    w <- as.numeric(w)
    w[!is.finite(w) | w <= 0] <- median(w[w > 0 & is.finite(w)], na.rm = TRUE)
    w[!is.finite(w) | w <= 0] <- 1
    sw <- sqrt(w)
    Aw <- A * sw
    bw <- b * sw
  }
  coef <- tryCatch(as.numeric(stats::lm.fit(Aw, bw)$coefficients),
                   error = function(e) rep(1 / ncol(A), ncol(A)))
  normalize_nonnegative(coef)
}

dwls_one_sample <- function(A, b, n.iter = 6, eps = 1e-8) {
  b <- as.numeric(b)
  b[!is.finite(b)] <- 0
  if (sum(b) > 0) b <- b / sum(b)
  x <- weighted_lm_nonnegative(A, b)
  for (it in seq_len(n.iter)) {
    fitted <- as.numeric(A %*% x)
    resid <- abs(b - fitted)
    scale <- stats::quantile(resid, 0.75, na.rm = TRUE) + eps
    w <- 1 / (resid + scale)
    cap <- stats::quantile(w, 0.95, na.rm = TRUE)
    w <- pmin(w, cap)
    w <- w / mean(w, na.rm = TRUE)
    x.new <- weighted_lm_nonnegative(A, b, w)
    if (max(abs(x.new - x)) < 1e-7) {
      x <- x.new
      break
    }
    x <- x.new
  }
  normalize_nonnegative(x)
}

run_dwls_matrix <- function(bulk, reference) {
  ref.prob <- normalize_cols(reference)
  bulk.prob <- normalize_cols(bulk)
  run_one <- function(j) dwls_one_sample(ref.prob, bulk.prob[, j])
  cl <- parallel::makeCluster(max(1, n_cores))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("ref.prob", "bulk.prob", "dwls_one_sample",
                                "weighted_lm_nonnegative", "normalize_nonnegative"),
                          envir = environment())
  fits <- parallel::parLapply(cl, seq_len(ncol(bulk.prob)), function(j) {
    dwls_one_sample(ref.prob, bulk.prob[, j])
  })
  mat <- do.call(rbind, fits)
  rownames(mat) <- colnames(bulk.prob)
  colnames(mat) <- colnames(ref.prob)
  mat
}

run_music_basic_matrix <- function(bulk, reference) {
  X <- normalize_cols(reference)
  Y <- normalize_cols(bulk)
  S <- rep(1, ncol(X))
  Sigma <- matrix(1e-8, nrow = nrow(X), ncol = ncol(X),
                  dimnames = list(rownames(X), colnames(X)))
  cl <- parallel::makeCluster(max(1, n_cores))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(MuSiC))
  parallel::clusterExport(cl, c("X", "Y", "S", "Sigma",
                                "weighted_lm_nonnegative", "normalize_nonnegative"),
                          envir = environment())
  fits <- parallel::parLapply(cl, seq_len(ncol(Y)), function(j) {
    fit <- tryCatch(
      MuSiC::music.basic(Y = as.numeric(Y[, j]), X = X, S = S, Sigma = Sigma,
                         iter.max = 1000, nu = 1e-4, eps = 0.01),
      error = function(e) NULL
    )
    if (is.null(fit) || any(!is.finite(fit$p.weight))) {
      weighted_lm_nonnegative(X, Y[, j])
    } else {
      normalize_nonnegative(fit$p.weight)
    }
  })
  mat <- do.call(rbind, fits)
  rownames(mat) <- colnames(Y)
  colnames(mat) <- colnames(X)
  mat
}

select_signature_genes <- function(reference, n = 1000) {
  ref <- log2(as.matrix(reference) + 1)
  state_var <- apply(ref, 1, stats::var, na.rm = TRUE)
  state_max <- apply(ref, 1, max, na.rm = TRUE)
  state_mean <- rowMeans(ref, na.rm = TRUE)
  score <- state_var * (state_max + 1e-8) / (state_mean + 1e-8)
  score[!is.finite(score)] <- 0
  genes <- names(sort(score, decreasing = TRUE))
  genes[seq_len(min(n, length(genes)))]
}

run_cibersort_nusvr_matrix <- function(bulk, reference) {
  genes <- select_signature_genes(reference, cibersort_signature_genes)
  X <- as.matrix(reference[genes, , drop = FALSE])
  Y <- as.matrix(bulk[genes, , drop = FALSE])
  X <- X[order(rownames(X)), , drop = FALSE]
  Y <- Y[order(rownames(Y)), , drop = FALSE]
  if (max(Y, na.rm = TRUE) < 50) Y <- 2^Y
  tmpc <- colnames(Y)
  tmpr <- rownames(Y)
  Y <- limma::normalizeQuantiles(Y)
  colnames(Y) <- tmpc
  rownames(Y) <- tmpr
  X <- (X - mean(X, na.rm = TRUE)) / stats::sd(as.vector(X), na.rm = TRUE)
  cibersort_core_one <- function(X, y) {
    nus <- c(0.25, 0.5, 0.75)
    fits <- lapply(nus, function(nu) {
      model <- e1071::svm(X, y, type = "nu-regression", kernel = "linear",
                          nu = nu, scale = FALSE)
      weights <- as.numeric(t(model$coefs) %*% model$SV)
      weights[weights < 0 | !is.finite(weights)] <- 0
      w <- normalize_nonnegative(weights)
      fitted <- as.numeric(X %*% w)
      list(w = w, rmse = sqrt(mean((fitted - y)^2, na.rm = TRUE)))
    })
    fits[[which.min(vapply(fits, function(z) z$rmse, numeric(1)))]]$w
  }
  cl <- parallel::makeCluster(max(1, n_cores))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(e1071))
  parallel::clusterExport(cl, c("X", "Y", "normalize_nonnegative", "cibersort_core_one"),
                          envir = environment())
  fits <- parallel::parLapply(cl, seq_len(ncol(Y)), function(j) {
    y <- Y[, j]
    y <- (y - mean(y, na.rm = TRUE)) / stats::sd(y, na.rm = TRUE)
    cibersort_core_one(X, y)
  })
  mat <- do.call(rbind, fits)
  rownames(mat) <- colnames(Y)
  colnames(mat) <- colnames(X)
  mat
}

run_citmic_native_matrix <- function(bulk) {
  citmic.input <- log2(as.matrix(bulk) + 1)
  pred <- CITMIC::CITMIC(as.data.frame(citmic.input), weighted = TRUE,
                         base = 10, damping = 0.90, cl.cores = n_cores)
  pred <- as.matrix(pred)
  if (!all(colnames(bulk) %in% rownames(pred)) && all(colnames(bulk) %in% colnames(pred))) {
    pred <- t(pred)
  }
  pred <- pred[colnames(bulk), , drop = FALSE]
  normalize_rows(pred)
}

run_immucellai_native_matrix <- function(bulk) {
  input_bulk <- bulk
  single_sample <- ncol(input_bulk) == 1
  original_sample <- colnames(input_bulk)
  if (single_sample) {
    input_bulk <- cbind(input_bulk, input_bulk[, 1, drop = FALSE])
    colnames(input_bulk) <- c(original_sample, paste0(original_sample, "_dup_for_ImmuCellAI"))
  }
  fit <- ImmuCellAI::ImmuCellAI_new(sample = as.data.frame(input_bulk),
                                   data_type = "RNA-seq",
                                   group_tag = FALSE,
                                   response_tag = FALSE,
                                   customer = 0)
  pred <- as.matrix(fit$Sample_abundance)
  pred <- pred[colnames(input_bulk), , drop = FALSE]
  if (single_sample) {
    pred <- pred[original_sample, , drop = FALSE]
  }
  normalize_rows(pred)
}

make_synthetic_bulk <- function(reference, n.samples, seed = 1) {
  set.seed(seed)
  n.states <- ncol(reference)
  alpha <- rep(0.7, n.states)
  prop <- matrix(stats::rgamma(n.samples * n.states, shape = rep(alpha, each = n.samples), rate = 1),
                 nrow = n.samples, ncol = n.states)
  prop <- prop / rowSums(prop)
  colnames(prop) <- colnames(reference)
  rownames(prop) <- sprintf("RuntimeSample_%03d", seq_len(n.samples))
  bulk <- as.matrix(reference) %*% t(prop)
  bulk <- bulk * 1e5 / pmax(colSums(bulk), 1e-12)
  colnames(bulk) <- rownames(prop)
  list(bulk = bulk, proportions = prop)
}

benchmark_method <- function(method, fun, genes.used, notes = "") {
  message("Running ", method, " ...")
  ensure_r_tempdir()
  gc()
  start <- Sys.time()
  cpu0 <- proc.time()
  status <- "OK"
  err <- ""
  pred <- tryCatch(fun(), error = function(e) {
    status <<- "FAILED"
    err <<- conditionMessage(e)
    NULL
  })
  cpu1 <- proc.time() - cpu0
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  if (!is.null(pred)) {
    ensure_r_tempdir()
    write.table(pred, file.path(out_dir, paste0(method, "_runtime_fraction.txt")),
                sep = "\t", quote = FALSE, row.names = TRUE)
  }
  data.frame(
    Method = method,
    Status = status,
    ElapsedSeconds = elapsed,
    UserCPUSeconds = as.numeric(cpu1["user.self"]),
    SystemCPUSeconds = as.numeric(cpu1["sys.self"]),
    SecondsPerSample = elapsed / n_benchmark_samples,
    Samples = n_benchmark_samples,
    GenesUsed = genes.used,
    CellStatesReturned = if (is.null(pred)) NA_integer_ else ncol(as.matrix(pred)),
    Cores = n_cores,
    Notes = notes,
    Error = err,
    stringsAsFactors = FALSE
  )
}

patch_immucellai_gsva()

message("Reading reference matrix...")
reference_all <- read_expression_matrix(reference_file)
hierarchy_all <- create_default_53_hierarchy(colnames(reference_all), hierarchy.mode = "tcell")
hierarchy_all <- ImmuCellAI2.0:::validate_hierarchy(hierarchy_all, reference_all)
reference_all <- reference_all[, hierarchy_all$state_name, drop = FALSE]

synthetic <- make_synthetic_bulk(reference_all, n_benchmark_samples, seed = seed)
bulk_all <- synthetic$bulk
write.table(synthetic$proportions, file.path(out_dir, "synthetic_true_proportions.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(bulk_all, file.path(out_dir, "synthetic_bulk_runtime_input.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

immune_genes <- intersect(read_gene_list(immune_gene_file), rownames(reference_all))
if (length(immune_genes) >= 100) {
  reference_immucellai2 <- reference_all[immune_genes, , drop = FALSE]
  bulk_immucellai2 <- bulk_all[immune_genes, , drop = FALSE]
  immucellai2_note <- "ImmuCellAI2 uses MarkerUsedDeconvolution immune genes"
} else {
  reference_immucellai2 <- reference_all
  bulk_immucellai2 <- bulk_all
  immucellai2_note <- "ImmuCellAI2 immune gene file unavailable; used all genes"
}

results <- list()

results[["ImmuCellAI2_tcell_VB"]] <- benchmark_method(
  "ImmuCellAI2_tcell_VB",
  function() {
    hierarchy <- create_default_53_hierarchy(colnames(reference_immucellai2), hierarchy.mode = "tcell")
    hierarchy <- ImmuCellAI2.0:::validate_hierarchy(hierarchy, reference_immucellai2)
    fit <- deconvolve_bulk_matrix(
      bulk.mat = bulk_immucellai2,
      reference = reference_immucellai2,
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
      seed = 123
    )
    fit$state.fraction
  },
  genes.used = nrow(bulk_immucellai2),
  notes = immucellai2_note
)

results[["DWLS_weighted_lm"]] <- benchmark_method(
  "DWLS_weighted_lm",
  function() run_dwls_matrix(bulk_all, reference_all),
  genes.used = nrow(bulk_all),
  notes = "weighted linear model implementation with 53-cell reference"
)

results[["MuSiC_basic"]] <- benchmark_method(
  "MuSiC_basic",
  function() run_music_basic_matrix(bulk_all, reference_all),
  genes.used = nrow(bulk_all),
  notes = "MuSiC::music.basic with 53-cell reference"
)

results[["CITMIC_native"]] <- benchmark_method(
  "CITMIC_native",
  function() run_citmic_native_matrix(bulk_all),
  genes.used = nrow(bulk_all),
  notes = "CITMIC native signatures; package API does not accept custom 53-cell reference"
)

results[["CIBERSORT_nuSVR_ref_top1000"]] <- benchmark_method(
  "CIBERSORT_nuSVR_ref_top1000",
  function() run_cibersort_nusvr_matrix(bulk_all, reference_all),
  genes.used = min(cibersort_signature_genes, nrow(bulk_all)),
  notes = paste0("nu-SVR on top ", cibersort_signature_genes, " discriminative reference genes")
)

results[["ImmuCellAI_native"]] <- benchmark_method(
  "ImmuCellAI_native",
  function() run_immucellai_native_matrix(bulk_all),
  genes.used = nrow(bulk_all),
  notes = "ImmuCellAI native signatures; package API does not accept custom 53-cell reference"
)

results[["BayesPrism_first_state_chain600_burn500"]] <- benchmark_method(
  "BayesPrism_first_state_chain600_burn500",
  function() {
    hierarchy_bp <- create_default_53_hierarchy(colnames(reference_all), hierarchy.mode = "tcell")
    hierarchy_bp <- ImmuCellAI2.0:::validate_hierarchy(hierarchy_bp, reference_all)
    reference.bp <- t(reference_all[, hierarchy_bp$state_name, drop = FALSE])
    mixture.bp <- t(bulk_all)
    set.seed(123)
    prism <- BayesPrism::new.prism(
      reference = reference.bp,
      input.type = "count.matrix",
      cell.type.labels = hierarchy_bp$major_lineage,
      cell.state.labels = hierarchy_bp$state_name,
      key = NA_character_,
      mixture = mixture.bp,
      outlier.cut = 0.01,
      outlier.fraction = 0.1,
      pseudo.min = 1e-8
    )
    bp <- fig3_run_prism(
      prism,
      n.cores = n_cores,
      update.gibbs = FALSE,
      gibbs.control = list(
        chain.length = chain_length,
        burn.in = burn_in,
        thinning = 2,
        alpha = 1,
        seed = 123
      )
    )
    BayesPrism::get.fraction(bp, which.theta = "first", state.or.type = "state")
  },
  genes.used = nrow(bulk_all),
  notes = paste0("chain.length=", chain_length, "; burn.in=", burn_in)
)

runtime_summary <- do.call(rbind, results)
runtime_summary <- runtime_summary[order(runtime_summary$ElapsedSeconds), ]
rownames(runtime_summary) <- NULL
write.table(runtime_summary, file.path(out_dir, "runtime_benchmark_7tools_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(out_dir, "runtime_benchmark_sessionInfo.txt"))

ensure_r_tempdir()
plot_df <- runtime_summary[runtime_summary$Status == "OK", , drop = FALSE]
plot_df$Method <- factor(plot_df$Method, levels = rev(plot_df$Method))
p <- ggplot(plot_df, aes(x = Method, y = ElapsedSeconds, fill = Method)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.35, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f s", ElapsedSeconds)), hjust = -0.08, size = 3.6) +
  coord_flip() +
  scale_fill_manual(values = rep(c("#F66463", "#379DA5", "#6BA5C3", "#FAC74C", "#95A8AC", "#AE997E", "#7B6DA8"), 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Elapsed seconds",
       title = paste0("Runtime benchmark, n=", n_benchmark_samples, ", cores=", n_cores)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.7),
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "runtime_benchmark_7tools_seconds_bar.png"), p, width = 10, height = 5, dpi = 300)
ggsave(file.path(out_dir, "runtime_benchmark_7tools_seconds_bar.pdf"), p, width = 10, height = 5)

print(runtime_summary)
