prepare_custom_expression <- function(x, label) {
  if (is.character(x) && length(x) == 1L) x <- read_expression_matrix(x)
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || any(dim(x) == 0L)) {
    stop("`", label, "` must be a non-empty numeric matrix or an expression file.")
  }
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(trimws(rownames(x)))) || any(!nzchar(trimws(colnames(x))))) {
    stop("`", label, "` must have non-empty gene row names and column names.")
  }
  if (anyDuplicated(colnames(x))) stop("`", label, "` column names must be unique.")
  if (any(!is.finite(x))) {
    warning("Non-finite values in `", label, "` were replaced by zero.")
    x[!is.finite(x)] <- 0
  }
  if (any(x < 0)) stop("`", label, "` must contain non-negative expression values.")
  if (anyDuplicated(rownames(x))) {
    warning("Duplicate gene identifiers in `", label, "` were collapsed by arithmetic mean.")
    x <- collapse_duplicate_genes(x, method = "mean")
  }
  x
}

#' Deconvolve with a user-supplied reference
#'
#' Uses all shared, informative genes and the deterministic flat Bayesian
#' posterior-mean engine. Reference columns are treated as separate populations
#' under a symmetric Dirichlet prior. No marker panel, biological hierarchy, or unknown
#' component is used. Arbitrary reference column names are supported.
#'
#' @param bulk Genes-by-samples non-negative numeric matrix or tab-delimited file.
#' @param reference Genes-by-populations non-negative numeric matrix or
#'   tab-delimited file. Gene identifiers must match those in bulk. TPM on a
#'   compatible scale is recommended for both inputs.
#' @param n.cores Number of samples processed concurrently.
#' @param pseudo.depth Pseudo-count depth after within-sample scaling.
#' @param n.iter Maximum number of deterministic posterior-mean updates.
#' @param vb.tol Convergence tolerance for the maximum fraction change.
#' @param alpha Total concentration of the symmetric Dirichlet prior.
#' @param seed Seed for worker initialization; the inference is deterministic.
#' @param verbose Print inference progress.
#' @return A list with subtype.fraction (samples by supplied reference columns),
#'   genes.used, inference.method, phi.final.list, and settings. No type or
#'   lineage aggregations are produced. Zero-expression samples and reference
#'   columns on retained genes are rejected rather than silently omitted.
#' @export
run_immucellai2_custom <- function(
  bulk,
  reference,
  n.cores = 1,
  pseudo.depth = 1e5,
  n.iter = 50,
  vb.tol = 1e-4,
  alpha = 1,
  seed = 123,
  verbose = FALSE
) {
  controls <- list(n.cores = n.cores, pseudo.depth = pseudo.depth,
                   n.iter = n.iter, vb.tol = vb.tol, alpha = alpha, seed = seed)
  for (name in names(controls)) {
    value <- controls[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("`", name, "` must be a finite numeric scalar.")
    }
  }
  if (n.cores < 1 || n.cores != floor(n.cores)) stop("`n.cores` must be a positive integer.")
  if (n.iter < 1 || n.iter != floor(n.iter)) stop("`n.iter` must be a positive integer.")
  if (pseudo.depth < 1) stop("`pseudo.depth` must be at least 1.")
  if (vb.tol < 0) stop("`vb.tol` must be non-negative.")
  if (alpha <= 0) stop("`alpha` must be positive.")
  if (seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    stop("`seed` must be a non-negative integer within R's integer range.")
  }

  bulk <- prepare_custom_expression(bulk, "bulk")
  reference <- prepare_custom_expression(reference, "reference")
  common <- intersect(rownames(bulk), rownames(reference))
  if (!length(common)) stop("No genes overlap bulk and reference. Check gene identifiers.")
  bulk <- bulk[common, , drop = FALSE]
  reference <- reference[common, , drop = FALSE]
  keep <- rowSums(bulk) > 0 & rowSums(reference) > 0
  if (!any(keep)) stop("No shared genes have positive expression in both inputs.")
  bulk <- bulk[keep, , drop = FALSE]
  reference <- reference[keep, , drop = FALSE]
  empty.samples <- colnames(bulk)[colSums(bulk) <= 0]
  empty.profiles <- colnames(reference)[colSums(reference) <= 0]
  if (length(empty.samples)) {
    stop("No expression on retained genes for bulk sample(s): ", paste(empty.samples, collapse = ", "))
  }
  if (length(empty.profiles)) {
    stop("No expression on retained genes for reference column(s): ", paste(empty.profiles, collapse = ", "))
  }
  count.mass <- vapply(seq_len(ncol(bulk)), function(j) {
    sum(round(bulk[, j] / sum(bulk[, j]) * pseudo.depth))
  }, numeric(1))
  if (any(count.mass <= 0)) stop("Pseudo-count rounding removed all evidence; increase `pseudo.depth`.")

  # Identity labels satisfy the matrix engine's schema without defining a hierarchy.
  populations <- colnames(reference)
  labels <- data.frame(state_name = populations, subtype = populations,
                       major_lineage = populations, stringsAsFactors = FALSE)
  fit <- deconvolve_bulk_matrix(
    bulk.mat = bulk, reference = reference, hierarchy = labels,
    hierarchy.mode = "flat", inference.method = "vb",
    add.unknown = FALSE, use.unknown.update = FALSE,
    pseudo.depth = pseudo.depth, n.iter = n.iter, vb.tol = vb.tol,
    alpha.major = alpha, n.cores = n.cores, seed = seed, verbose = verbose
  )
  result <- fit[c("subtype.fraction", "genes.used", "inference.method", "phi.final.list")]
  result$settings <- list(
    workflow = "custom_reference", inference.method = "vb",
    gene.selection = "all_shared_informative", genes.shared = length(common),
    genes.used = length(fit$genes.used), reference.populations = populations,
    n.cores = n.cores, pseudo.depth = pseudo.depth, n.iter = n.iter,
    vb.tol = vb.tol, alpha = alpha, seed = seed
  )
  result
}
