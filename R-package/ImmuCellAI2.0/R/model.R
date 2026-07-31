make_initial_theta_generic <- function(major.lineage.labels, subtype.labels, state.labels) {
  lineages <- unique(major.lineage.labels)
  theta.major <- rep(1 / length(lineages), length(lineages))
  names(theta.major) <- lineages

  theta.sub <- list()
  for (lin in lineages) {
    idx.lin <- which(major.lineage.labels == lin)
    subtype.names <- unique(subtype.labels[idx.lin])
    theta.sub[[lin]] <- rep(1 / length(subtype.names), length(subtype.names))
    names(theta.sub[[lin]]) <- subtype.names
  }

  theta.state <- list()
  for (sbt in unique(subtype.labels)) {
    idx.sbt <- which(subtype.labels == sbt)
    state.names <- state.labels[idx.sbt]
    theta.state[[sbt]] <- rep(1 / length(state.names), length(state.names))
    names(theta.state[[sbt]]) <- state.names
  }

  list(theta.major = theta.major, theta.sub = theta.sub, theta.state = theta.state)
}

update_theta_hierarchical <- function(Z, major.lineage.labels, subtype.labels, state.labels,
                                      alpha.major = 10, alpha.sub = 5, alpha.state = 1) {
  if (is.null(colnames(Z))) colnames(Z) <- state.labels
  lineages <- unique(major.lineage.labels)
  subtypes <- unique(subtype.labels)
  state.counts <- colSums(Z)
  state.counts <- state.counts[state.labels]

  lineage.counts <- sapply(lineages, function(lin) {
    idx <- which(major.lineage.labels == lin)
    sum(state.counts[idx])
  })
  theta.major <- sample_dirichlet_safe(alpha.major + lineage.counts)
  names(theta.major) <- lineages

  theta.sub <- list()
  for (lin in lineages) {
    idx.lin <- which(major.lineage.labels == lin)
    subtype.names <- unique(subtype.labels[idx.lin])
    subtype.counts <- sapply(subtype.names, function(sbt) {
      idx <- which(subtype.labels == sbt)
      sum(state.counts[idx])
    })
    theta.sub[[lin]] <- sample_dirichlet_safe(alpha.sub + subtype.counts)
    names(theta.sub[[lin]]) <- subtype.names
  }

  theta.state <- list()
  for (sbt in subtypes) {
    idx.sbt <- which(subtype.labels == sbt)
    state.names <- state.labels[idx.sbt]
    local.counts <- state.counts[idx.sbt]
    theta.state[[sbt]] <- sample_dirichlet_safe(alpha.state + local.counts)
    names(theta.state[[sbt]]) <- state.names
  }

  theta.final <- rep(0, length(state.labels))
  names(theta.final) <- state.labels
  for (i in seq_along(state.labels)) {
    st <- state.labels[i]
    sbt <- subtype.labels[i]
    lin <- major.lineage.labels[i]
    theta.final[st] <- theta.major[lin] * theta.sub[[lin]][sbt] * theta.state[[sbt]][st]
  }
  theta.final <- normalize_prob(theta.final)

  list(theta.major = theta.major, theta.sub = theta.sub, theta.state = theta.state,
       theta.final = theta.final)
}

update_Z_hierarchical <- function(bulk.vec, phi, theta.major, theta.sub, theta.state,
                                  major.lineage.labels, subtype.labels, state.labels) {
  Z <- matrix(0, nrow = nrow(phi), ncol = ncol(phi))
  rownames(Z) <- rownames(phi)
  colnames(Z) <- colnames(phi)

  state.prior <- rep(0, length(state.labels))
  names(state.prior) <- state.labels
  for (i in seq_along(state.labels)) {
    st <- state.labels[i]
    sbt <- subtype.labels[i]
    lin <- major.lineage.labels[i]
    state.prior[st] <- theta.major[lin] * theta.sub[[lin]][sbt] * theta.state[[sbt]][st]
  }
  state.prior <- normalize_prob(state.prior)

  for (g in seq_len(nrow(phi))) {
    y_g <- bulk.vec[g]
    if (is.na(y_g) || y_g <= 0) next
    prob <- phi[g, ] * state.prior[colnames(phi)]
    prob <- normalize_prob(prob)
    Z[g, ] <- as.numeric(stats::rmultinom(n = 1, size = y_g, prob = prob))
  }
  Z
}

make_state_prior_from_hierarchy <- function(major.lineage.labels, subtype.labels, state.labels) {
  prior <- rep(0, length(state.labels))
  names(prior) <- state.labels
  lineages <- unique(major.lineage.labels)
  major.prior <- rep(1 / length(lineages), length(lineages))
  names(major.prior) <- lineages
  for (lin in lineages) {
    idx.lin <- which(major.lineage.labels == lin)
    subtypes <- unique(subtype.labels[idx.lin])
    sub.prior <- rep(1 / length(subtypes), length(subtypes))
    names(sub.prior) <- subtypes
    for (sbt in subtypes) {
      idx.sbt <- idx.lin[subtype.labels[idx.lin] == sbt]
      prior[idx.sbt] <- major.prior[lin] * sub.prior[sbt] / length(idx.sbt)
    }
  }
  normalize_prob(prior)
}

posterior_mean_theta_hierarchical <- function(
  state.counts,
  major.lineage.labels,
  subtype.labels,
  state.labels,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  prior.theta = NULL
) {
  state.counts <- state.counts[state.labels]
  state.counts[is.na(state.counts)] <- 0

  if (is.null(prior.theta)) {
    prior.theta <- make_state_prior_from_hierarchy(major.lineage.labels, subtype.labels, state.labels)
  } else {
    prior.theta <- normalize_prob(prior.theta[state.labels])
    names(prior.theta) <- state.labels
  }

  lineages <- unique(major.lineage.labels)
  subtypes <- unique(subtype.labels)

  lineage.counts <- sapply(lineages, function(lin) {
    sum(state.counts[major.lineage.labels == lin])
  })
  prior.major <- sapply(lineages, function(lin) {
    sum(prior.theta[major.lineage.labels == lin])
  })
  prior.major <- normalize_prob(prior.major)
  theta.major <- normalize_prob(alpha.major * prior.major + lineage.counts)
  names(theta.major) <- lineages

  theta.sub <- list()
  for (lin in lineages) {
    idx.lin <- which(major.lineage.labels == lin)
    subtype.names <- unique(subtype.labels[idx.lin])
    subtype.counts <- sapply(subtype.names, function(sbt) {
      sum(state.counts[subtype.labels == sbt])
    })
    prior.sub <- sapply(subtype.names, function(sbt) {
      sum(prior.theta[subtype.labels == sbt])
    })
    prior.sub <- normalize_prob(prior.sub)
    theta.sub[[lin]] <- normalize_prob(alpha.sub * prior.sub + subtype.counts)
    names(theta.sub[[lin]]) <- subtype.names
  }

  theta.state <- list()
  for (sbt in subtypes) {
    idx.sbt <- which(subtype.labels == sbt)
    state.names <- state.labels[idx.sbt]
    local.counts <- state.counts[idx.sbt]
    prior.state <- normalize_prob(prior.theta[idx.sbt])
    theta.state[[sbt]] <- normalize_prob(alpha.state * prior.state + local.counts)
    names(theta.state[[sbt]]) <- state.names
  }

  theta.final <- rep(0, length(state.labels))
  names(theta.final) <- state.labels
  for (i in seq_along(state.labels)) {
    st <- state.labels[i]
    sbt <- subtype.labels[i]
    lin <- major.lineage.labels[i]
    theta.final[st] <- theta.major[lin] * theta.sub[[lin]][sbt] * theta.state[[sbt]][st]
  }
  theta.final <- normalize_prob(theta.final)

  list(
    theta.major = theta.major,
    theta.sub = theta.sub,
    theta.state = theta.state,
    theta.final = theta.final
  )
}

expected_state_counts <- function(bulk.vec, phi, theta.final) {
  theta.final <- normalize_prob(theta.final[colnames(phi)])
  weighted <- sweep(phi, 2, theta.final, "*")
  denom <- rowSums(weighted)
  denom[is.na(denom) | denom <= 0] <- 1e-300
  resp <- weighted / denom
  bulk.vec <- as.numeric(bulk.vec[rownames(phi)])
  bulk.vec[is.na(bulk.vec) | bulk.vec < 0] <- 0
  counts <- colSums(resp * bulk.vec)
  names(counts) <- colnames(phi)
  list(counts = counts, loglik = sum(bulk.vec * log(denom)))
}

#' Run vectorized hierarchical variational Bayes for one bulk sample
#'
#' This keeps the original major -> subtype -> state Dirichlet hierarchy, but
#' replaces stochastic multinomial allocation with expected latent counts. It is
#' therefore deterministic, faster than Gibbs sampling, and still returns a
#' hierarchical Bayesian posterior mean for state fractions.
#'
#' @param bulk.vec Named numeric vector of bulk expression/counts.
#' @param phi.ref Gene x state reference probability matrix.
#' @param hierarchy data.frame with state_name, subtype, major_lineage.
#' @param unknown.state Name of UNKNOWN state.
#' @param n.iter Maximum variational iterations.
#' @param burn.in Delay before optional UNKNOWN-profile adaptation; ordinary VB
#'   posterior means are not averaged after burn-in.
#' @param tol Stop when max absolute theta change is below this value.
#' @param alpha.major,alpha.sub,alpha.state Hierarchical Dirichlet prior strengths.
#' @param prior.theta Optional state-level prior probabilities.
#' @param use.unknown.update Adapt the UNKNOWN profile from residual expression.
#' @param adaptive.params List overriding residual-update settings.
#' @param verbose Print iteration diagnostics.
#' @return list with theta.posterior, theta.trace, loglik.trace, and phi.final.
#' @export
run_hierarchical_vb_model <- function(
  bulk.vec,
  phi.ref,
  hierarchy,
  unknown.state = "UNKNOWN",
  n.iter = 80,
  burn.in = 0,
  tol = 1e-7,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  prior.theta = NULL,
  use.unknown.update = FALSE,
  adaptive.params = list(),
  verbose = FALSE
) {
  state.labels <- hierarchy$state_name
  subtype.labels <- hierarchy$subtype
  major.lineage.labels <- hierarchy$major_lineage

  phi.current <- phi.ref[, state.labels, drop = FALSE]
  bulk.vec <- bulk.vec[rownames(phi.current)]
  bulk.vec[is.na(bulk.vec)] <- 0
  bulk.vec[bulk.vec < 0] <- 0

  if (is.null(prior.theta)) {
    prior.theta <- make_state_prior_from_hierarchy(major.lineage.labels, subtype.labels, state.labels)
  } else {
    prior.theta <- normalize_prob(prior.theta[state.labels])
    names(prior.theta) <- state.labels
  }

  theta.final <- prior.theta
  theta.trace <- matrix(0, nrow = n.iter, ncol = length(state.labels))
  colnames(theta.trace) <- state.labels
  loglik.trace <- rep(NA_real_, n.iter)

  default.adaptive <- default_adaptive_params()
  adaptive <- modifyList(default.adaptive, adaptive.params)

  last.iter <- n.iter
  for (iter in seq_len(n.iter)) {
    ez <- expected_state_counts(bulk.vec, phi.current, theta.final)
    theta.obj <- posterior_mean_theta_hierarchical(
      state.counts = ez$counts,
      major.lineage.labels = major.lineage.labels,
      subtype.labels = subtype.labels,
      state.labels = state.labels,
      alpha.major = alpha.major,
      alpha.sub = alpha.sub,
      alpha.state = alpha.state,
      prior.theta = prior.theta
    )

    delta <- max(abs(theta.obj$theta.final[state.labels] - theta.final[state.labels]))
    theta.final <- theta.obj$theta.final[state.labels]
    theta.trace[iter, ] <- theta.final
    loglik.trace[iter] <- ez$loglik

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      cat("VB iteration:", iter, "delta:", signif(delta, 4), "loglik:", signif(ez$loglik, 6), "\n")
    }

    if (use.unknown.update && unknown.state %in% state.labels && iter > burn.in) {
      phi.current <- do.call(
        update_phi_unknown_residual_gated_fallback,
        c(list(
          bulk.vec = bulk.vec,
          phi.current = phi.current,
          theta.final = theta.final,
          unknown.state = unknown.state
        ), adaptive)
      )
    }

    if (delta < tol) {
      last.iter <- iter
      break
    }
  }

  theta.trace <- theta.trace[seq_len(last.iter), , drop = FALSE]
  loglik.trace <- loglik.trace[seq_len(last.iter)]
  list(
    theta.posterior = theta.final,
    theta.trace = theta.trace,
    loglik.trace = loglik.trace,
    phi.final = phi.current,
    n.iter.used = last.iter
  )
}

#' Run flat variational Bayes for one bulk sample
#'
#' This is a non-hierarchical Bayesian deconvolution model. It uses a single
#' Dirichlet prior over all states and expected latent counts for deterministic
#' posterior updates. It is useful when the hierarchy is too conservative or
#' when the reference itself generated the mixtures.
#'
#' @param bulk.vec Named numeric vector of bulk expression/counts.
#' @param phi.ref Gene x state reference probability matrix.
#' @param hierarchy data.frame with state_name, subtype, major_lineage.
#' @param unknown.state Name of the optional residual state.
#' @param n.iter Maximum variational iterations.
#' @param burn.in Delay before optional UNKNOWN-profile adaptation; ordinary VB
#'   posterior means are not averaged after burn-in.
#' @param tol Stop when max absolute theta change is below this value.
#' @param alpha.state Dirichlet prior strength over states.
#' @param prior.theta Optional state-level prior probabilities.
#' @param use.unknown.update Adapt the UNKNOWN profile from residual expression.
#' @param adaptive.params List overriding residual-update settings.
#' @param verbose Print iteration diagnostics.
#' @return list with theta.posterior, theta.trace, loglik.trace, and phi.final.
#' @export
run_flat_vb_model <- function(
  bulk.vec,
  phi.ref,
  hierarchy,
  unknown.state = "UNKNOWN",
  n.iter = 80,
  burn.in = 0,
  tol = 1e-7,
  alpha.state = 1,
  prior.theta = NULL,
  use.unknown.update = FALSE,
  adaptive.params = list(),
  verbose = FALSE
) {
  state.labels <- hierarchy$state_name
  phi.current <- phi.ref[, state.labels, drop = FALSE]
  bulk.vec <- bulk.vec[rownames(phi.current)]
  bulk.vec[is.na(bulk.vec)] <- 0
  bulk.vec[bulk.vec < 0] <- 0

  if (is.null(prior.theta)) {
    prior.theta <- rep(1 / length(state.labels), length(state.labels))
    names(prior.theta) <- state.labels
  } else {
    prior.theta <- normalize_prob(prior.theta[state.labels])
    names(prior.theta) <- state.labels
  }

  theta.final <- prior.theta
  theta.trace <- matrix(0, nrow = n.iter, ncol = length(state.labels))
  colnames(theta.trace) <- state.labels
  loglik.trace <- rep(NA_real_, n.iter)
  adaptive <- modifyList(default_adaptive_params(), adaptive.params)

  last.iter <- n.iter
  for (iter in seq_len(n.iter)) {
    ez <- expected_state_counts(bulk.vec, phi.current, theta.final)
    theta.new <- normalize_prob(alpha.state * prior.theta + ez$counts[state.labels])
    names(theta.new) <- state.labels

    delta <- max(abs(theta.new - theta.final))
    theta.final <- theta.new
    theta.trace[iter, ] <- theta.final
    loglik.trace[iter] <- ez$loglik

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      cat("Flat VB iteration:", iter, "delta:", signif(delta, 4), "loglik:", signif(ez$loglik, 6), "\n")
    }

    if (use.unknown.update && unknown.state %in% state.labels && iter > burn.in) {
      phi.current <- do.call(
        update_phi_unknown_residual_gated_fallback,
        c(list(
          bulk.vec = bulk.vec,
          phi.current = phi.current,
          theta.final = theta.final,
          unknown.state = unknown.state
        ), adaptive)
      )
    }

    if (delta < tol) {
      last.iter <- iter
      break
    }
  }

  theta.trace <- theta.trace[seq_len(last.iter), , drop = FALSE]
  loglik.trace <- loglik.trace[seq_len(last.iter)]
  list(
    theta.posterior = theta.final,
    theta.trace = theta.trace,
    loglik.trace = loglik.trace,
    phi.final = phi.current,
    n.iter.used = last.iter
  )
}

update_phi_unknown_residual_gated_fallback <- function(
  bulk.vec,
  phi.current,
  theta.final,
  unknown.state = "UNKNOWN",
  lambda.ref = 5000,
  learning.rate = 0.05,
  residual.quantile = 0.50,
  min.unknown.theta = 0.03,
  min.residual.fraction = 0.02,
  max.effective.residual.genes = 150,
  strong.residual.fraction = 0.05,
  strong.max.effective.genes = 120,
  fallback.learning.rate = 0.01,
  fallback.lambda.ref = 10000
) {
  phi.updated <- phi.current
  if (!unknown.state %in% colnames(phi.current)) return(phi.updated)

  known.states <- setdiff(colnames(phi.current), unknown.state)
  known.theta <- theta.final[known.states]
  known.theta[is.na(known.theta)] <- 0

  known.mean <- as.numeric(phi.current[, known.states, drop = FALSE] %*% known.theta)
  known.mean <- known.mean * sum(bulk.vec)

  residual <- bulk.vec - known.mean
  residual[is.na(residual)] <- 0
  residual[residual < 0] <- 0
  positive.residual <- residual[residual > 0]
  if (length(positive.residual) == 0) return(phi.updated)

  cutoff <- stats::quantile(positive.residual, probs = residual.quantile, na.rm = TRUE)
  residual[residual < cutoff] <- 0
  if (sum(residual) <= 0) return(phi.updated)

  residual.fraction <- sum(residual) / sum(bulk.vec)
  residual.profile <- residual / sum(residual)
  residual.entropy <- -sum(residual.profile * log(residual.profile + 1e-12))
  residual.effective.genes <- exp(residual.entropy)

  theta.unknown <- theta.final[unknown.state]
  if (is.na(theta.unknown)) theta.unknown <- 0

  unknown.theta.ok <- theta.unknown >= min.unknown.theta
  residual.evidence.ok <- residual.fraction >= min.residual.fraction &&
    residual.effective.genes <= max.effective.residual.genes
  strong.residual.evidence.ok <- residual.fraction >= strong.residual.fraction &&
    residual.effective.genes <= strong.max.effective.genes

  if (unknown.theta.ok && residual.evidence.ok) {
    effective.lambda <- lambda.ref
    effective.learning.rate <- learning.rate
  } else if (strong.residual.evidence.ok) {
    effective.lambda <- fallback.lambda.ref
    effective.learning.rate <- fallback.learning.rate
  } else {
    return(phi.updated)
  }

  phi.posterior <- effective.lambda * phi.current[, unknown.state] + residual
  phi.posterior <- phi.posterior / sum(phi.posterior)
  phi.updated[, unknown.state] <- (1 - effective.learning.rate) * phi.current[, unknown.state] +
    effective.learning.rate * phi.posterior
  phi.updated[, unknown.state] <- phi.updated[, unknown.state] / sum(phi.updated[, unknown.state])
  phi.updated
}

#' Run gated-fallback adaptive hierarchical model for one bulk sample
#'
#' @param bulk.vec Named numeric vector of bulk expression/counts.
#' @param phi.ref Gene x state reference probability matrix. Columns should sum to 1.
#' @param hierarchy data.frame with state_name, subtype, major_lineage.
#' @param unknown.state Name of UNKNOWN state.
#' @param n.iter Number of legacy MCMC iterations.
#' @param burn.in Number of initial iterations discarded from posterior averaging.
#' @param alpha.major,alpha.sub,alpha.state Hierarchical Dirichlet prior strengths.
#' @param use.unknown.update Whether to update UNKNOWN using gated-fallback residual adaptation.
#' @param adaptive.params List overriding residual-update settings.
#' @param verbose Print iteration diagnostics.
#' @return list with theta.posterior, theta.trace, and phi.final.
#' @export
run_gated_fallback_adaptive_model <- function(
  bulk.vec,
  phi.ref,
  hierarchy,
  unknown.state = "UNKNOWN",
  n.iter = 150,
  burn.in = 50,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  use.unknown.update = TRUE,
  adaptive.params = list(),
  verbose = FALSE
) {
  state.labels <- hierarchy$state_name
  subtype.labels <- hierarchy$subtype
  major.lineage.labels <- hierarchy$major_lineage

  phi.ref <- phi.ref[, state.labels, drop = FALSE]
  bulk.vec <- bulk.vec[rownames(phi.ref)]
  bulk.vec[is.na(bulk.vec)] <- 0
  bulk.vec[bulk.vec < 0] <- 0

  init <- make_initial_theta_generic(major.lineage.labels, subtype.labels, state.labels)
  theta.major <- init$theta.major
  theta.sub <- init$theta.sub
  theta.state <- init$theta.state
  phi.current <- phi.ref

  theta.trace <- matrix(0, nrow = n.iter, ncol = length(state.labels))
  colnames(theta.trace) <- state.labels

  default.adaptive <- list(
    lambda.ref = 5000,
    learning.rate = 0.05,
    residual.quantile = 0.50,
    min.unknown.theta = 0.03,
    min.residual.fraction = 0.02,
    max.effective.residual.genes = 150,
    strong.residual.fraction = 0.05,
    strong.max.effective.genes = 120,
    fallback.learning.rate = 0.01,
    fallback.lambda.ref = 10000
  )
  adaptive <- modifyList(default.adaptive, adaptive.params)

  for (iter in seq_len(n.iter)) {
    if (verbose && iter %% 25 == 0) cat("Iteration:", iter, "\n")

    Z <- update_Z_hierarchical(
      bulk.vec = bulk.vec,
      phi = phi.current,
      theta.major = theta.major,
      theta.sub = theta.sub,
      theta.state = theta.state,
      major.lineage.labels = major.lineage.labels,
      subtype.labels = subtype.labels,
      state.labels = state.labels
    )

    theta.obj <- update_theta_hierarchical(
      Z = Z,
      major.lineage.labels = major.lineage.labels,
      subtype.labels = subtype.labels,
      state.labels = state.labels,
      alpha.major = alpha.major,
      alpha.sub = alpha.sub,
      alpha.state = alpha.state
    )

    theta.major <- theta.obj$theta.major
    theta.sub <- theta.obj$theta.sub
    theta.state <- theta.obj$theta.state
    theta.trace[iter, ] <- theta.obj$theta.final[state.labels]

    if (use.unknown.update && unknown.state %in% state.labels && iter > burn.in) {
      phi.current <- do.call(
        update_phi_unknown_residual_gated_fallback,
        c(list(
          bulk.vec = bulk.vec,
          phi.current = phi.current,
          theta.final = theta.obj$theta.final,
          unknown.state = unknown.state
        ), adaptive)
      )
    }
  }

  theta.posterior <- colMeans(theta.trace[(burn.in + 1):n.iter, , drop = FALSE])
  list(theta.posterior = theta.posterior, theta.trace = theta.trace, phi.final = phi.current)
}
