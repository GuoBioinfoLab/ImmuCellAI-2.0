#' Normalize matrix columns to sum to one
#'
#' @param mat Numeric matrix.
#' @return Column-normalized matrix.
normalize_columns_fast <- function(mat) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  s <- colSums(mat)
  s[is.na(s) | s <= 0] <- 1
  sweep(mat, 2, s, "/")
}

#' Select state-specific marker genes from a reference matrix
#'
#' Markers are selected by log1p expression above each gene's median reference
#' expression across states. This is intentionally simple and deterministic so
#' benchmark sweeps are reproducible.
#'
#' @param reference Gene x state reference expression matrix.
#' @param marker.n Number of markers per state.
#' @return Named list of marker gene vectors.
#' @export
select_state_markers <- function(reference, marker.n = 50) {
  reference <- as.matrix(reference)
  mode(reference) <- "numeric"
  reference[is.na(reference)] <- 0
  reference[reference < 0] <- 0

  lr <- log1p(reference)
  med <- stats::median
  gene.median <- apply(lr, 1, med, na.rm = TRUE)
  score <- sweep(lr, 1, gene.median, "-")
  marker.n <- max(1, as.integer(marker.n))

  out <- vector("list", ncol(reference))
  names(out) <- colnames(reference)
  for (j in seq_len(ncol(reference))) {
    ord <- order(score[, j], decreasing = TRUE, na.last = NA)
    out[[j]] <- rownames(reference)[ord[seq_len(min(marker.n, length(ord)))]]
  }
  out
}

#' Build marker-score abundance estimates
#'
#' This estimator is optimized for increasing-proportion spike-in simulations:
#' each state's abundance signal is the mean normalized bulk expression of that
#' state's marker genes. The raw score is useful for monotonicity/correlation
#' tests; optional column normalization returns fraction-like outputs.
#'
#' @param bulk.mat Gene x sample bulk expression matrix.
#' @param reference Gene x state reference expression matrix.
#' @param hierarchy Optional hierarchy table.
#' @param marker.n Number of markers per state.
#' @param normalize.scores If TRUE, normalize state scores within each sample.
#' @return Deconvolution-like result list with state.score and state.fraction.
#' @export
deconvolve_marker_score <- function(
  bulk.mat,
  reference,
  hierarchy = NULL,
  marker.n = 50,
  normalize.scores = TRUE
) {
  bulk.mat <- collapse_duplicate_genes(as.matrix(bulk.mat), method = "mean")
  reference <- collapse_duplicate_genes(as.matrix(reference), method = "mean")
  mode(bulk.mat) <- "numeric"
  mode(reference) <- "numeric"
  bulk.mat[is.na(bulk.mat)] <- 0
  reference[is.na(reference)] <- 0
  bulk.mat[bulk.mat < 0] <- 0
  reference[reference < 0] <- 0

  common.genes <- intersect(rownames(bulk.mat), rownames(reference))
  bulk <- bulk.mat[common.genes, , drop = FALSE]
  ref <- reference[common.genes, , drop = FALSE]
  keep <- rowSums(bulk) > 0 & rowSums(ref) > 0
  bulk <- bulk[keep, , drop = FALSE]
  ref <- ref[keep, , drop = FALSE]

  if (is.null(hierarchy)) {
    hierarchy <- create_default_53_hierarchy(colnames(ref), hierarchy.mode = "hybrid")
  }
  hierarchy <- validate_hierarchy(hierarchy, ref)
  ref <- ref[, hierarchy$state_name, drop = FALSE]

  bulk.norm <- normalize_columns_fast(bulk)
  markers <- select_state_markers(ref, marker.n = marker.n)
  score <- matrix(0, nrow = ncol(bulk.norm), ncol = ncol(ref))
  rownames(score) <- colnames(bulk.norm)
  colnames(score) <- colnames(ref)

  for (state in colnames(ref)) {
    genes <- intersect(markers[[state]], rownames(bulk.norm))
    if (length(genes) > 0) {
      score[, state] <- colMeans(bulk.norm[genes, , drop = FALSE])
    }
  }

  fraction <- score
  if (isTRUE(normalize.scores)) {
    rs <- rowSums(fraction)
    rs[is.na(rs) | rs <= 0] <- 1
    fraction <- sweep(fraction, 1, rs, "/")
  }

  ag <- aggregate_by_hierarchy(fraction, hierarchy)
  list(
    state.fraction = fraction,
    state.score = score,
    subtype.fraction = ag$subtype,
    major.fraction = ag$major,
    hierarchy = hierarchy,
    marker.n = marker.n,
    genes.used = rownames(bulk),
    markers = markers
  )
}

project_simplex <- function(v) {
  v <- as.numeric(v)
  if (length(v) == 1) return(1)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  ind <- seq_along(v)
  keep <- u - cssv / ind > 0
  if (!any(keep)) return(rep(1 / length(v), length(v)))
  theta <- cssv[max(which(keep))] / max(which(keep))
  w <- pmax(v - theta, 0)
  s <- sum(w)
  if (s <= 0) rep(1 / length(v), length(v)) else w / s
}

solve_simplex_pgd <- function(ref, bulk, n.iter = 80) {
  ref <- as.matrix(ref)
  bulk <- as.matrix(bulk)
  gram <- crossprod(ref)
  rhs <- crossprod(ref, bulk)
  eig <- suppressWarnings(eigen(gram, symmetric = TRUE, only.values = TRUE)$values)
  step <- 1 / max(max(eig, na.rm = TRUE), 1e-12)
  w <- matrix(1 / ncol(ref), nrow = ncol(ref), ncol = ncol(bulk))
  for (iter in seq_len(n.iter)) {
    w <- w - step * (gram %*% w - rhs)
    for (j in seq_len(ncol(w))) w[, j] <- project_simplex(w[, j])
  }
  rownames(w) <- colnames(ref)
  colnames(w) <- colnames(bulk)
  w
}

union_marker_genes <- function(reference, marker.n) {
  unique(unlist(select_state_markers(reference, marker.n = marker.n), use.names = FALSE))
}

aggregate_reference_by_label <- function(reference, labels) {
  labels <- as.character(labels)
  groups <- unique(labels)
  out <- matrix(0, nrow = nrow(reference), ncol = length(groups))
  rownames(out) <- rownames(reference)
  colnames(out) <- groups
  for (g in groups) {
    out[, g] <- rowMeans(reference[, labels == g, drop = FALSE])
  }
  out
}

fit_flat_simplex <- function(reference, bulk, marker.n = 400, n.iter = 80) {
  genes <- union_marker_genes(reference, marker.n)
  genes <- intersect(genes, rownames(reference))
  ref.x <- normalize_columns_fast(reference[genes, , drop = FALSE])
  bulk.x <- normalize_columns_fast(bulk[genes, , drop = FALSE])
  solve_simplex_pgd(ref.x, bulk.x, n.iter = n.iter)
}

#' Fast deterministic hierarchical deconvolution
#'
#' Uses marker selection plus projected-gradient simplex fitting. In
#' hierarchical mode, major lineages, subtypes, and states are fitted in a
#' top-down product, which is much faster and more reproducible than MCMC for
#' parameter sweeps.
#'
#' @param bulk.mat Gene x sample bulk expression matrix.
#' @param reference Gene x state reference expression matrix.
#' @param hierarchy Optional hierarchy table.
#' @param hierarchy.mode One of flat, hierarchical, hybrid, or tcell. The legacy value tcell_only is accepted as an alias.
#' @param marker.n Number of markers per state.
#' @param n.iter Projected-gradient iterations.
#' @return Deconvolution result list.
#' @export
deconvolve_bulk_fast <- function(
  bulk.mat,
  reference,
  hierarchy = NULL,
  hierarchy.mode = c("hybrid", "hierarchical", "flat", "tcell", "tcell_only"),
  marker.n = 400,
  n.iter = 80
) {
  hierarchy.mode <- normalize_hierarchy_mode(hierarchy.mode)
  bulk.mat <- collapse_duplicate_genes(as.matrix(bulk.mat), method = "mean")
  reference <- collapse_duplicate_genes(as.matrix(reference), method = "mean")
  common.genes <- intersect(rownames(bulk.mat), rownames(reference))
  bulk <- bulk.mat[common.genes, , drop = FALSE]
  ref <- reference[common.genes, , drop = FALSE]
  mode(bulk) <- "numeric"
  mode(ref) <- "numeric"
  bulk[is.na(bulk)] <- 0
  ref[is.na(ref)] <- 0
  bulk[bulk < 0] <- 0
  ref[ref < 0] <- 0
  keep <- rowSums(bulk) > 0 & rowSums(ref) > 0
  bulk <- bulk[keep, , drop = FALSE]
  ref <- ref[keep, , drop = FALSE]

  if (is.null(hierarchy)) {
    hierarchy <- create_default_53_hierarchy(colnames(ref), hierarchy.mode = hierarchy.mode)
  }
  hierarchy <- validate_hierarchy(hierarchy, ref)
  states <- hierarchy$state_name
  ref <- ref[, states, drop = FALSE]

  if (hierarchy.mode == "flat") {
    w <- t(fit_flat_simplex(ref, bulk, marker.n = marker.n, n.iter = n.iter))
    colnames(w) <- states
    rownames(w) <- colnames(bulk)
  } else {
    major.ref <- aggregate_reference_by_label(ref, hierarchy$major_lineage)
    major.w <- fit_flat_simplex(major.ref, bulk, marker.n = marker.n, n.iter = n.iter)
    final <- matrix(0, nrow = length(states), ncol = ncol(bulk))
    rownames(final) <- states
    colnames(final) <- colnames(bulk)

    for (major in colnames(major.ref)) {
      idx.major <- which(hierarchy$major_lineage == major)
      sub.ref <- aggregate_reference_by_label(ref[, idx.major, drop = FALSE], hierarchy$subtype[idx.major])
      sub.w <- fit_flat_simplex(sub.ref, bulk, marker.n = marker.n, n.iter = n.iter)
      for (subtype in colnames(sub.ref)) {
        idx.sub <- idx.major[hierarchy$subtype[idx.major] == subtype]
        if (length(idx.sub) == 1) {
          local.w <- matrix(1, nrow = 1, ncol = ncol(bulk))
          rownames(local.w) <- states[idx.sub]
          colnames(local.w) <- colnames(bulk)
        } else {
          local.w <- fit_flat_simplex(ref[, idx.sub, drop = FALSE], bulk, marker.n = marker.n, n.iter = n.iter)
        }
        final[states[idx.sub], ] <- sweep(local.w, 2, major.w[major, ] * sub.w[subtype, ], "*")
      }
    }

    cs <- colSums(final)
    cs[is.na(cs) | cs <= 0] <- 1
    w <- t(sweep(final, 2, cs, "/"))
  }

  ag <- aggregate_by_hierarchy(w, hierarchy)
  list(
    state.fraction = w,
    subtype.fraction = ag$subtype,
    major.fraction = ag$major,
    hierarchy = hierarchy,
    hierarchy.mode = hierarchy.mode,
    marker.n = marker.n,
    genes.used = rownames(bulk)
  )
}
