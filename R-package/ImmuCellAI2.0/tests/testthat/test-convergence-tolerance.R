convergence_fixture <- function() {
  set.seed(173)
  ref <- matrix(rexp(120 * 4, rate = 10), nrow = 120)
  for (j in seq_len(4)) ref[(30 * (j - 1) + 1):(30 * j), j] <-
    ref[(30 * (j - 1) + 1):(30 * j), j] + 5
  dimnames(ref) <- list(paste0("Gene", seq_len(120)), paste0("Cell", seq_len(4)))
  phi <- sweep(ref, 2, colSums(ref), "/")
  counts <- as.numeric(round(1e5 * phi %*% c(0.55, 0.25, 0.15, 0.05)))
  names(counts) <- rownames(ref)
  hierarchy <- data.frame(state_name = colnames(ref), subtype = c("A", "A", "B", "B"),
    major_lineage = c("L1", "L1", "L2", "L2"), stringsAsFactors = FALSE)
  list(ref = ref, phi = phi, counts = counts, hierarchy = hierarchy)
}

test_that("all VB entry points default to the same tolerance", {
  expect_equal(formals(deconvolve_bulk_matrix)$vb.tol, 1e-4)
  expect_equal(formals(run_hierarchical_vb_model)$tol, 1e-4)
  expect_equal(formals(run_flat_vb_model)$tol, 1e-4)
  expect_equal(formals(run_immucellai2)$vb.tol, 1e-4)
  expect_equal(immucellai2_defaults()$vb.tol, 1e-4)
  expect_equal(formals(run_immucellai2)$n.iter, 50)
  expect_equal(immucellai2_defaults()$n.iter, 50)
})

test_that("VB stops below the new tolerance and still honors overrides and caps", {
  d <- convergence_fixture()
  for (runner in list(run_hierarchical_vb_model, run_flat_vb_model)) {
    arguments <- list(bulk.vec = d$counts, phi.ref = d$phi, hierarchy = d$hierarchy, n.iter = 500)
    default <- do.call(runner, arguments)
    explicit <- do.call(runner, c(arguments, list(tol = 1e-4)))
    strict <- do.call(runner, c(arguments, list(tol = 1e-6)))
    expect_equal(default$theta.posterior, explicit$theta.posterior, tolerance = 0)
    expect_equal(default$n.iter.used, explicit$n.iter.used)
    expect_gt(default$n.iter.used, 1)
    expect_lt(default$n.iter.used, 500)
    delta <- apply(abs(diff(default$theta.trace)), 1, max)
    expect_lt(tail(delta, 1), 1e-4)
    expect_true(all(head(delta, -1) >= 1e-4))
    expect_lte(default$n.iter.used, strict$n.iter.used)
    expect_lt(max(abs(diff(tail(strict$theta.trace, 2)))), 1e-6)
    capped <- do.call(runner, list(bulk.vec = d$counts, phi.ref = d$phi,
      hierarchy = d$hierarchy, n.iter = 1, tol = 0))
    expect_equal(capped$n.iter.used, 1)
  }
})

test_that("matrix and standard wrappers propagate the new default", {
  d <- convergence_fixture()
  bulk <- matrix(d$counts, ncol = 1, dimnames = list(names(d$counts), "S1"))
  for (method in c("vb", "flat_vb")) {
    arguments <- list(bulk.mat = bulk, reference = d$ref, hierarchy = d$hierarchy,
      inference.method = method, n.iter = 50, n.cores = 1)
    default <- suppressWarnings(do.call(deconvolve_bulk_matrix, arguments))
    explicit <- suppressWarnings(do.call(deconvolve_bulk_matrix, c(arguments, list(vb.tol = 1e-4))))
    expect_equal(default$subtype.fraction, explicit$subtype.fraction, tolerance = 0)
  }
  default <- suppressWarnings(run_immucellai2(bulk, reference = d$ref,
    marker.genes = rownames(d$ref), hierarchy.mode = "flat"))
  strict <- suppressWarnings(run_immucellai2(bulk, reference = d$ref,
    marker.genes = rownames(d$ref), hierarchy.mode = "flat", vb.tol = 1e-6))
  expect_equal(default$settings$vb.tol, 1e-4)
  expect_equal(default$settings$n.iter, 50)
  expect_equal(strict$settings$vb.tol, 1e-6)
})
