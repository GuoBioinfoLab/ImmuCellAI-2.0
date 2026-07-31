add_unknown_to_reference <- function(reference, hierarchy, unknown.state = "UNKNOWN", unknown.mode = c("uniform", "mean")) {
  unknown.mode <- match.arg(unknown.mode)
  if (unknown.state %in% colnames(reference)) {
    if (!unknown.state %in% hierarchy$state_name) {
      hierarchy <- rbind(hierarchy, data.frame(
        state_name = unknown.state,
        subtype = "UNKNOWN",
        major_lineage = "UNKNOWN",
        stringsAsFactors = FALSE
      ))
    }
    return(list(reference = reference, hierarchy = hierarchy))
  }
  if (unknown.mode == "uniform") {
    unknown.profile <- rep(1, nrow(reference))
  } else {
    unknown.profile <- rowMeans(reference)
  }
  reference <- cbind(reference, unknown.profile)
  colnames(reference)[ncol(reference)] <- unknown.state
  hierarchy <- rbind(hierarchy, data.frame(
    state_name = unknown.state,
    subtype = "UNKNOWN",
    major_lineage = "UNKNOWN",
    stringsAsFactors = FALSE
  ))
  list(reference = reference, hierarchy = hierarchy)
}

#' Deconvolve a bulk expression matrix
#'
#' Lower-level engine supporting deterministic hierarchical or flat VB and a
#' legacy MCMC compatibility path. New analyses should normally use
#' run_immucellai2().
#'
#' @param bulk.mat Gene-by-sample non-negative expression matrix.
#' @param reference Gene-by-state non-negative reference matrix.
#' @param hierarchy Optional table with state_name, subtype, and major_lineage.
#' @param hierarchy.mode One of flat, hierarchical, hybrid, or tcell_only.
#' @param inference.method One of vb, flat_vb, or the legacy mcmc path.
#' @param add.unknown Add an optional latent residual state.
#' @param unknown.state Name assigned to the residual state.
#' @param unknown.mode Initial residual profile: uniform or reference mean.
#' @param pseudo.depth Total pseudo-count mass used for each sample.
#' @param n.iter Maximum number of VB updates or legacy MCMC iterations.
#' @param burn.in Used by the legacy MCMC path and optional residual adaptation;
#'   it does not define posterior averaging for ordinary VB runs.
#' @param alpha.major,alpha.sub,alpha.state Dirichlet prior concentrations.
#' @param vb.tol Maximum absolute fraction change used as the VB stopping rule.
#' @param prior.theta Optional named state-level prior probability vector.
#' @param use.unknown.update Adapt the UNKNOWN profile from residual expression.
#' @param adaptive.params Optional list overriding residual-update parameters.
#' @param n.cores Number of samples processed concurrently.
#' @param seed Reproducibility seed used for workers and legacy sampling.
#' @param verbose Print per-sample or per-iteration progress.
#' @return A list containing sample-by-state, subtype, and major fraction
#'   matrices, the hierarchy, genes used, and per-sample final references.
#' @export
 deconvolve_bulk_matrix <- function(
  bulk.mat,
  reference,
 hierarchy = NULL,
  hierarchy.mode = c("flat", "hierarchical", "hybrid", "tcell_only"),
  inference.method = c("mcmc", "vb", "flat_vb"),
  add.unknown = FALSE,
  unknown.state = "UNKNOWN",
  unknown.mode = c("uniform", "mean"),
  pseudo.depth = 1e5,
  n.iter = 150,
  burn.in = 50,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  vb.tol = 1e-7,
  prior.theta = NULL,
  use.unknown.update = add.unknown,
  adaptive.params = list(),
  n.cores = 1,
  seed = 123,
  verbose = FALSE
) {
  unknown.mode <- match.arg(unknown.mode)
  hierarchy.mode <- match.arg(hierarchy.mode)
  inference.method <- match.arg(inference.method)
  set.seed(seed)

  bulk.mat <- as.matrix(bulk.mat)
  reference <- as.matrix(reference)
  mode(bulk.mat) <- "numeric"
  mode(reference) <- "numeric"
  bulk.mat[is.na(bulk.mat)] <- 0
  bulk.mat[bulk.mat < 0] <- 0
  reference[is.na(reference)] <- 0
  reference[reference < 0] <- 0

  bulk.mat <- collapse_duplicate_genes(bulk.mat, method = "mean")
  reference <- collapse_duplicate_genes(reference, method = "mean")

  if (is.null(hierarchy)) hierarchy <- create_default_53_hierarchy(colnames(reference), hierarchy.mode = hierarchy.mode)
  hierarchy <- validate_hierarchy(hierarchy, reference)

  if (add.unknown) {
    tmp <- add_unknown_to_reference(reference, hierarchy, unknown.state, unknown.mode)
    reference <- tmp$reference
    hierarchy <- tmp$hierarchy
  }

  state.labels <- hierarchy$state_name
  common.genes <- intersect(rownames(bulk.mat), rownames(reference))
  if (length(common.genes) < 1000) warning("Common genes fewer than 1000: ", length(common.genes))

  bulk <- bulk.mat[common.genes, , drop = FALSE]
  ref <- reference[common.genes, state.labels, drop = FALSE]
  keep <- rowSums(bulk) > 0 & rowSums(ref) > 0
  bulk <- bulk[keep, , drop = FALSE]
  ref <- ref[keep, , drop = FALSE]

  phi.ref <- safe_col_normalize(ref + 1e-8)

  run_one <- function(sample.name) {
    bulk.vec <- bulk[, sample.name]
    bulk.vec[is.na(bulk.vec)] <- 0
    bulk.vec[bulk.vec < 0] <- 0
    if (sum(bulk.vec) <= 0) return(NULL)
    bulk.count <- round(bulk.vec / sum(bulk.vec) * pseudo.depth)
    names(bulk.count) <- rownames(bulk)
    if (inference.method == "flat_vb") {
      fit <- run_flat_vb_model(
        bulk.vec = bulk.count,
        phi.ref = phi.ref,
        hierarchy = hierarchy,
        unknown.state = unknown.state,
        n.iter = n.iter,
        burn.in = burn.in,
        tol = vb.tol,
        use.unknown.update = use.unknown.update,
        alpha.state = alpha.state,
        prior.theta = prior.theta,
        adaptive.params = adaptive.params,
        verbose = verbose
      )
    } else if (inference.method == "vb") {
      fit <- run_hierarchical_vb_model(
        bulk.vec = bulk.count,
        phi.ref = phi.ref,
        hierarchy = hierarchy,
        unknown.state = unknown.state,
        n.iter = n.iter,
        burn.in = burn.in,
        tol = vb.tol,
        use.unknown.update = use.unknown.update,
        alpha.major = alpha.major,
        alpha.sub = alpha.sub,
        alpha.state = alpha.state,
        prior.theta = prior.theta,
        adaptive.params = adaptive.params,
        verbose = verbose
      )
    } else {
      fit <- run_gated_fallback_adaptive_model(
        bulk.vec = bulk.count,
        phi.ref = phi.ref,
        hierarchy = hierarchy,
        unknown.state = unknown.state,
        n.iter = n.iter,
        burn.in = burn.in,
        use.unknown.update = use.unknown.update,
        alpha.major = alpha.major,
        alpha.sub = alpha.sub,
        alpha.state = alpha.state,
        adaptive.params = adaptive.params,
        verbose = verbose
      )
    }
    list(sample = sample.name, theta = fit$theta.posterior[state.labels], phi.final = fit$phi.final)
  }

  sample.names <- colnames(bulk)

  if (n.cores > 1) {
    n.cores.use <- min(n.cores, parallel::detectCores() - 1)
    n.cores.use <- max(1, n.cores.use)
    cl <- parallel::makeCluster(n.cores.use)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c(
        "bulk", "phi.ref", "hierarchy", "state.labels", "unknown.state", "pseudo.depth",
        "n.iter", "burn.in", "alpha.major", "alpha.sub", "alpha.state", "vb.tol", "prior.theta",
        "inference.method", "use.unknown.update", "adaptive.params", "verbose",
        "normalize_prob", "safe_col_normalize", "sample_dirichlet_safe",
        "make_initial_theta_generic", "update_theta_hierarchical", "update_Z_hierarchical",
        "make_state_prior_from_hierarchy", "posterior_mean_theta_hierarchical", "expected_state_counts",
        "update_phi_unknown_residual_gated_fallback", "run_gated_fallback_adaptive_model",
        "default_adaptive_params", "run_hierarchical_vb_model", "run_flat_vb_model", "run_one"
      ),
      envir = environment()
    )
    parallel::clusterSetRNGStream(cl, seed)
    fits <- parallel::parLapply(cl, sample.names, run_one)
  } else {
    fits <- lapply(sample.names, function(smp) {
      if (verbose) cat("Running sample:", smp, "\n")
      run_one(smp)
    })
  }

  fits <- fits[!sapply(fits, is.null)]
  theta.list <- lapply(fits, function(x) x$theta)
  names(theta.list) <- sapply(fits, function(x) x$sample)
  state.mat <- do.call(rbind, theta.list)

  ag <- aggregate_by_hierarchy(state.mat, hierarchy)

  list(
    state.fraction = state.mat,
    subtype.fraction = ag$subtype,
    major.fraction = ag$major,
    hierarchy = hierarchy,
    hierarchy.mode = hierarchy.mode,
    inference.method = inference.method,
    genes.used = rownames(bulk),
    phi.final.list = lapply(fits, function(x) x$phi.final)
  )
}
