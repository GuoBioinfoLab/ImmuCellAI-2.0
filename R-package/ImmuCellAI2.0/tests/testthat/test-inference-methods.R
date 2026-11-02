inference_method_fixture <- function() {
  set.seed(93)
  ref <- matrix(rexp(1200 * 3), nrow = 1200,
                dimnames = list(paste0("Gene", 1:1200), c("CellA", "CellB", "CellC")))
  bulk <- ref %*% cbind(S1 = c(.6, .3, .1), S2 = c(.2, .3, .5))
  list(ref = ref, bulk = bulk)
}

test_that("the matrix API exposes only vb and mcmc", {
  expect_identical(eval(formals(deconvolve_bulk_matrix)$inference.method), c("mcmc", "vb"))
  expect_error(deconvolve_bulk_matrix(NULL, NULL, inference.method = "flat_vb"), "mcmc.*vb")
  expect_false("inference.method" %in% names(formals(run_immucellai2)))
  expect_identical(immucellai2_defaults()$inference.method, "vb")
})

test_that("automatic flat updates preserve the independent-branch posterior", {
  x <- inference_method_fixture()
  for (unknown in c(FALSE, TRUE)) {
    ref <- x$ref
    if (unknown) ref <- cbind(ref, UNKNOWN = 1)
    states <- colnames(ref)
    h <- data.frame(state_name = states, subtype = states, major_lineage = states)
    phi <- sweep(ref + 1e-8, 2, colSums(ref + 1e-8), "/")
    prior <- setNames(seq_along(states) / sum(seq_along(states)), states)
    for (alpha in c(.1, 1, 5)) {
      fit <- deconvolve_bulk_matrix(x$bulk, ref, hierarchy = h, hierarchy.mode = "flat",
        inference.method = "vb", add.unknown = unknown, alpha.major = alpha,
        prior.theta = prior, n.iter = 10, burn.in = 2)
      for (sample in colnames(x$bulk)) {
        counts <- round(x$bulk[, sample] / sum(x$bulk[, sample]) * 1e5)
        expected <- run_hierarchical_vb_model(counts, phi, h, alpha.major = alpha,
          prior.theta = prior, n.iter = 10, burn.in = 2, use.unknown.update = unknown)
        expect_equal(fit$subtype.fraction[sample, ], expected$theta.posterior,
                     tolerance = 1e-12)
        expect_equal(fit$phi.final.list[[match(sample, colnames(x$bulk))]],
                     expected$phi.final, tolerance = 1e-12)
      }
    }
  }
})

test_that("both hierarchy modes propagate to parallel VB workers", {
  x <- inference_method_fixture()
  colnames(x$ref) <- c("CD4Tcm", "CD4Tem", "BGC")
  for (mode in c("tcell", "flat")) {
    for (unknown in c(FALSE, TRUE)) {
      args <- list(bulk.mat = x$bulk, reference = x$ref, hierarchy.mode = mode,
                   inference.method = "vb", add.unknown = unknown, n.iter = 3)
      single <- do.call(deconvolve_bulk_matrix, c(args, list(n.cores = 1)))
      multi <- do.call(deconvolve_bulk_matrix, c(args, list(n.cores = 2)))
      expect_identical(multi, single)
    }
  }
})
