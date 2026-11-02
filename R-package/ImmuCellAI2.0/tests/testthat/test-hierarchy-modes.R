test_that("only tcell and flat are accepted at hierarchy entry points", {
  states <- c("CD4Tcm", "CD4Tem", "CD8Tcm", "BGC", "cMo")
  expect_identical(create_default_53_hierarchy(states, "tcell"),
                   create_tcell_53_hierarchy(states))
  expect_identical(create_default_53_hierarchy(states, "flat"),
                   create_flat_53_hierarchy(states))
  expect_equal(immucellai2_defaults()$hierarchy.mode, "tcell")

  for (mode in c("hierarchical", "hybrid", "tcell_only", "invalid")) {
    expect_error(create_default_53_hierarchy(states, mode), "flat.*tcell")
    expect_error(run_immucellai2(NULL, hierarchy.mode = mode), "flat.*tcell")
    expect_error(deconvolve_bulk_matrix(NULL, NULL, hierarchy.mode = mode), "flat.*tcell")
    expect_error(deconvolve_bulk_fast(NULL, NULL, hierarchy.mode = mode), "flat.*tcell")
  }
  expect_false("refine.tcell.subtypes" %in% names(formals(create_default_53_hierarchy)))
  expect_false("include.hybrid" %in% names(formals(run_experiment1_benchmark)))
  retired <- c("create_hybrid_53_hierarchy", "create_tcell_only_53_hierarchy",
               "refine_tcell_subtypes", "refine_monocyte_subtypes")
  expect_false(any(retired %in% getNamespaceExports("ImmuCellAI2.0")))
})

hierarchy_mode_fixture <- function() {
  set.seed(62)
  ref <- matrix(rexp(120 * 5), nrow = 120)
  rownames(ref) <- paste0("Gene", seq_len(nrow(ref)))
  colnames(ref) <- c("CD4Tcm", "CD4Tem", "CD4Tn", "BGC", "cMo")
  truth <- c(.1, .2, .3)
  weights <- rbind(truth, matrix(rep((1 - truth) / 4, each = 4), nrow = 4))
  bulk <- ref %*% weights
  colnames(bulk) <- paste0("CD4Tcm_", seq_along(truth))
  info <- data.frame(Sample = colnames(bulk), TargetCell = "CD4Tcm", Truth = truth)
  list(ref = ref, bulk = bulk, info = info)
}

test_that("auxiliary estimators default to the supported tcell mode", {
  x <- hierarchy_mode_fixture()
  h <- create_tcell_53_hierarchy(colnames(x$ref))
  expect_identical(deconvolve_marker_score(x$bulk, x$ref, marker.n = 10)$hierarchy, h)
  fit <- deconvolve_bulk_fast(x$bulk, x$ref, marker.n = 10, n.iter = 3)
  expect_identical(fit$hierarchy.mode, "tcell")
  expect_identical(fit$hierarchy, h)
  expect_identical(fit, deconvolve_bulk_fast(x$bulk, x$ref, hierarchy.mode = "tcell",
                                           marker.n = 10, n.iter = 3))
})

test_that("comparison helpers run only supported hierarchy configurations", {
  x <- hierarchy_mode_fixture()
  out <- tempfile("two-hierarchy-modes-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- suppressMessages(run_experiment1_benchmark(
    x$bulk, x$ref, out.dir = file.path(out, "benchmark"), sample.info = x$info,
    n.iter = 4, burn.in = 2))
  expect_setequal(result$methods$hierarchy.mode, c("flat", "tcell"))
  expect_setequal(names(result$results), c("flat_no_UNKNOWN", "flat_UNKNOWN",
                                          "tcell_no_UNKNOWN", "tcell_UNKNOWN"))
  expect_equal(nrow(result$methods), 4L)
  for (i in seq_len(nrow(result$methods))) {
    fit <- result$results[[result$methods$Method[i]]]
    expect_identical(fit$hierarchy.mode, result$methods$hierarchy.mode[i])
  }
  grid <- suppressMessages(run_fast_hierarchy_grid(
    x$bulk, x$ref, out.dir = file.path(out, "grid"), marker.grid = 10,
    n.iter = 3, sample.info = x$info))
  expect_setequal(names(grid$metrics), c("marker_score_raw_10", "marker_score_fraction_10",
                                       "simplex_flat_10", "simplex_tcell_10"))
})
