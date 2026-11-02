custom_reference_fixture <- function() {
  set.seed(93)
  ref <- matrix(rexp(1200 * 3), nrow = 1200)
  for (k in 1:3) ref[(400 * (k - 1) + 1):(400 * k), k] <-
    ref[(400 * (k - 1) + 1):(400 * k), k] + 12
  dimnames(ref) <- list(paste0("CustomGene", 1:1200), c("Population A", "Population-B", "Population.C"))
  truth <- rbind(S1 = c(.6, .3, .1), S2 = c(.2, .3, .5))
  bulk <- sweep(ref, 2, colSums(ref), "/") %*% t(truth) * 1e6
  list(ref = ref, bulk = bulk, truth = truth)
}

test_that("custom reference uses every shared valid gene and no bundled resources", {
  x <- custom_reference_fixture()
  local_mocked_bindings(
    load_immucellai2_reference = function(...) stop("bundled atlas called"),
    load_immucellai2_markers = function(...) stop("bundled markers called"),
    create_default_53_hierarchy = function(...) stop("default hierarchy called"),
    run_hierarchical_vb_model = function(...) stop("hierarchical inference called"),
    update_phi_unknown_residual_gated_fallback = function(...) stop("UNKNOWN update called")
  )
  ref <- rbind(x$ref, OnlyReference = c(4, 5, 6), ZeroReference = c(0, 0, 0),
               ZeroBulk = c(1, 2, 3))
  bulk <- rbind(x$bulk, OnlyBulk = c(2, 3), ZeroReference = c(4, 5), ZeroBulk = c(0, 0))
  fit <- run_immucellai2_custom(bulk, ref)
  expect_identical(fit$genes.used, rownames(x$ref))
  expect_identical(colnames(fit$subtype.fraction), colnames(ref))
  expect_identical(rownames(fit$subtype.fraction), colnames(bulk))
  expect_equal(dim(fit$subtype.fraction), c(2L, 3L))
  expect_true(all(is.finite(fit$subtype.fraction)))
  expect_true(all(fit$subtype.fraction >= 0))
  expect_equal(unname(rowSums(fit$subtype.fraction)), c(1, 1), tolerance = 1e-12)
  expect_lt(max(abs(unname(fit$subtype.fraction) - unname(x$truth))), 0.01)
  expect_false(any(c("hierarchy", "type.fraction", "lineage.fraction") %in% names(fit)))
  expect_identical(fit$inference.method, "vb")
  expect_identical(fit$settings$gene.selection, "all_shared_informative")
  expect_equal(fit$settings$genes.used, 1200L)
  expect_equal(fit$settings$genes.shared, 1202L)
  expect_identical(run_immucellai2(bulk, reference = ref), fit)
})

test_that("custom updates match the existing flat engine and posterior mean formula", {
  x <- custom_reference_fixture()
  labels <- data.frame(state_name = colnames(x$ref), subtype = colnames(x$ref),
                       major_lineage = colnames(x$ref))
  for (iterations in c(1L, 50L)) {
    fit <- run_immucellai2_custom(x$bulk, x$ref, n.iter = iterations)
    core <- deconvolve_bulk_matrix(x$bulk, x$ref, hierarchy = labels,
      hierarchy.mode = "flat", inference.method = "vb", add.unknown = FALSE,
      alpha.major = 1, use.unknown.update = FALSE, n.iter = iterations, vb.tol = 1e-4)
    expect_identical(fit$subtype.fraction, core$subtype.fraction)
    expect_identical(fit$phi.final.list, core$phi.final.list)
  }
  phi <- sweep(x$ref + 1e-8, 2, colSums(x$ref + 1e-8), "/")
  responsibility <- phi / rowSums(phi)
  counts <- round(x$bulk[, 1] / sum(x$bulk[, 1]) * 1e5)
  allocated <- colSums(responsibility * counts)
  expected <- (allocated + 1 / ncol(phi)) / (sum(allocated) + 1)
  one <- run_immucellai2_custom(x$bulk[, 1, drop = FALSE], x$ref, n.iter = 1)
  expect_equal(one$subtype.fraction[1, ], expected, tolerance = 1e-14)
  fit <- run_immucellai2_custom(x$bulk, x$ref)
  other.seed <- run_immucellai2_custom(x$bulk, x$ref, seed = 8)
  expect_identical(fit$subtype.fraction, other.seed$subtype.fraction)
  multi <- run_immucellai2_custom(x$bulk, x$ref, n.cores = 2)
  expect_identical(fit$subtype.fraction, multi$subtype.fraction)
})

test_that("custom reference does not expose marker, hierarchy or unknown options", {
  forbidden <- c("marker.genes", "hierarchy", "hierarchy.mode", "add.unknown", "unknown.mode")
  expect_false(any(forbidden %in% names(formals(run_immucellai2_custom))))
  x <- custom_reference_fixture()
  for (option in list(list(marker.genes = rownames(x$ref)), list(hierarchy.mode = "flat"),
                     list(add.unknown = TRUE), list(add.unknown = FALSE))) {
    expect_error(do.call(run_immucellai2, c(list(bulk = x$bulk, reference = x$ref), option)),
                 "Custom references use all shared genes")
    expect_error(do.call(run_immucellai2_custom, c(list(bulk = x$bulk, reference = x$ref), option)),
                 "unused argument")
  }
})

test_that("file inputs, arbitrary labels, gene order and custom export are supported", {
  x <- custom_reference_fixture()
  out <- tempfile("custom-input-")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  utils::write.table(x$bulk, file.path(out, "bulk.txt"), sep = "\t", quote = FALSE, col.names = NA)
  utils::write.table(x$ref, file.path(out, "reference.txt"), sep = "\t", quote = FALSE, col.names = NA)
  from.files <- run_immucellai2_custom(file.path(out, "bulk.txt"), file.path(out, "reference.txt"))
  fit <- run_immucellai2_custom(x$bulk, x$ref)
  expect_equal(from.files$subtype.fraction, fit$subtype.fraction, tolerance = 1e-12)
  reordered <- run_immucellai2_custom(x$bulk[1200:1, , drop = FALSE], x$ref[c(2:1200, 1), ])
  expect_equal(reordered$subtype.fraction, fit$subtype.fraction, tolerance = 1e-12)
  single <- run_immucellai2_custom(x$bulk[, 1, drop = FALSE], x$ref[, 1, drop = FALSE])
  expect_equal(dim(single$subtype.fraction), c(1L, 1L))
  expect_equal(unname(single$subtype.fraction[1, 1]), 1)
  write_deconvolution_outputs(fit, file.path(out, "output"))
  expect_identical(list.files(file.path(out, "output")), "subtype_fraction.txt")
})

test_that("custom inputs reject absent evidence and invalid dimensions or controls", {
  x <- custom_reference_fixture()
  invalid <- x$bulk
  invalid[, 1] <- 0
  expect_error(run_immucellai2_custom(invalid, x$ref), "bulk sample.*S1")
  invalid <- x$ref
  invalid[, 1] <- 0
  expect_error(run_immucellai2_custom(x$bulk, invalid), "reference column.*Population A")
  rownames(invalid) <- paste0("absent", seq_len(nrow(invalid)))
  expect_error(run_immucellai2_custom(x$bulk, invalid), "No genes overlap")
  expect_error(run_immucellai2_custom(x$bulk * 0, x$ref), "No shared genes")
  expect_error(run_immucellai2_custom(-x$bulk, x$ref), "non-negative")
  expect_error(run_immucellai2_custom(unname(x$bulk), x$ref), "non-empty gene row names")
  invalid <- x$ref
  colnames(invalid)[2] <- colnames(invalid)[1]
  expect_error(run_immucellai2_custom(x$bulk, invalid), "column names must be unique")
  expect_error(run_immucellai2_custom(x$bulk, x$ref, n.iter = 0), "positive integer")
  expect_error(run_immucellai2_custom(x$bulk, x$ref, n.cores = 1.5), "positive integer")
  expect_error(run_immucellai2_custom(x$bulk, x$ref, alpha = 0), "positive")
  expect_error(run_immucellai2_custom(x$bulk, x$ref, vb.tol = NA_real_), "finite numeric scalar")
  expect_error(run_immucellai2_custom(x$bulk, x$ref, pseudo.depth = 1), "rounding removed all evidence")
})
