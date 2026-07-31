normalize_prob <- function(x) {
  x[is.na(x)] <- 0
  x[x < 0] <- 0
  s <- sum(x)
  if (s <= 0) return(rep(1 / length(x), length(x)))
  x / s
}

safe_col_normalize <- function(mat) {
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  for (j in seq_len(ncol(mat))) {
    s <- sum(mat[, j])
    if (s <= 0) {
      mat[, j] <- rep(1 / nrow(mat), nrow(mat))
    } else {
      mat[, j] <- mat[, j] / s
    }
  }
  mat
}

sample_dirichlet_safe <- function(alpha.vec) {
  alpha.vec[is.na(alpha.vec)] <- 1e-8
  alpha.vec[alpha.vec < 1e-8] <- 1e-8
  x <- stats::rgamma(n = length(alpha.vec), shape = alpha.vec, rate = 1)
  x / sum(x)
}

collapse_duplicate_genes <- function(mat, method = c("mean", "sum")) {
  method <- match.arg(method)
  mat <- as.matrix(mat)
  gene <- rownames(mat)
  if (!any(duplicated(gene))) return(mat)
  mat2 <- rowsum(mat, group = gene)
  if (method == "mean") {
    counts <- table(gene)
    mat2 <- sweep(mat2, 1, as.numeric(counts[rownames(mat2)]), "/")
  }
  mat2
}
