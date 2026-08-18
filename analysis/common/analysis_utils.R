options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Install them before running this analysis."
    )
  }
}

assert_file <- function(path, label = basename(path)) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop(label, " was not found: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_table_auto <- function(path) {
  assert_file(path)
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    require_packages("readxl")
    return(as.data.frame(readxl::read_excel(path), check.names = FALSE))
  }
  require_packages("data.table")
  data.table::fread(path, data.table = FALSE, check.names = FALSE)
}

read_fraction_matrix <- function(path) {
  x <- read_table_auto(path)
  if (ncol(x) < 2L) stop("Fraction table must contain a sample column and cell fractions.")
  sample_id <- as.character(x[[1L]])
  x[[1L]] <- NULL
  mat <- as.matrix(x)
  storage.mode(mat) <- "double"
  rownames(mat) <- sample_id
  if (any(!is.finite(mat))) stop("Fraction table contains non-finite values.")
  mat
}

write_matrix_with_id <- function(mat, path, id_name = "sample") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- data.frame(sample = rownames(mat), mat, check.names = FALSE)
  names(out)[1L] <- id_name
  data.table::fwrite(out, path, sep = "\t", quote = FALSE)
}

run_standard_deconvolution <- function(expression_file, output_dir, n_cores = 1L) {
  require_packages(c("ImmuCellAI2.0", "data.table"))
  expression_file <- assert_file(expression_file, "Expression matrix")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- ImmuCellAI2.0::run_immucellai2(
    bulk = expression_file,
    hierarchy.mode = "tcell_only",
    add.unknown = FALSE,
    n.cores = as.integer(n_cores),
    pseudo.depth = 1e5,
    n.iter = 50,
    vb.tol = 1e-6,
    alpha.major = 10,
    alpha.sub = 5,
    alpha.state = 1,
    seed = 123
  )
  ImmuCellAI2.0::write_deconvolution_outputs(fit, output_dir)
  saveRDS(fit, file.path(output_dir, "immucellai2_result.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(fit)
}

zscore_rows <- function(sample_by_feature, clip = 1) {
  z <- t(scale(sample_by_feature))
  z[!is.finite(z)] <- 0
  z[z > clip] <- clip
  z[z < -clip] <- -clip
  z
}

clr_transform <- function(x, eps = 1e-5) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  lx <- log(pmax(x, 0) + eps)
  sweep(lx, 1, rowMeans(lx), "-")
}

auc_rank <- function(y, score) {
  ok <- !is.na(y) & is.finite(score)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  if (length(unique(y)) != 2L) return(NA_real_)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  (sum(rank(score)[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_points <- function(y, score) {
  ok <- !is.na(y) & is.finite(score)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  data.frame(
    FPR = c(0, cumsum(y == 0L) / sum(y == 0L), 1),
    TPR = c(0, cumsum(y == 1L) / sum(y == 1L), 1)
  )
}

save_plot_pair <- function(plot, output_stem, width, height, dpi = 600) {
  ggplot2::ggsave(paste0(output_stem, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(output_stem, ".png"), plot, width = width, height = height, dpi = dpi)
}
