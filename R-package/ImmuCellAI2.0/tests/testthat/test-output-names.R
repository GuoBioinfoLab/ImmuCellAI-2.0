output_fixture <- function() {
  set.seed(42)
  ref <- matrix(rexp(1000 * 5), nrow = 1000)
  rownames(ref) <- paste0("Gene", seq_len(nrow(ref)))
  colnames(ref) <- c("CD4Tcm", "CD4Tem", "CD4Tn", "CD8Tcm", "BGC")
  bulk <- ref %*% cbind(S1 = c(.1, .2, .3, .1, .3), S2 = c(.3, .1, .1, .2, .3))
  list(bulk = bulk, ref = ref)
}

expect_fraction_levels <- function(fit) {
  expect_true(all(c("subtype.fraction", "type.fraction", "lineage.fraction") %in% names(fit)))
  expect_false(any(c("state.fraction", "major.fraction") %in% names(fit)))
  expect_false(anyDuplicated(names(fit)) > 0L)
  h <- fit$hierarchy
  terminal <- fit$subtype.fraction
  expect_identical(colnames(terminal), h$state_name)
  for (field in c("subtype.fraction", "type.fraction", "lineage.fraction")) {
    mat <- fit[[field]]
    expect_true(all(is.finite(mat)))
    expect_true(all(mat >= 0))
    expect_identical(rownames(mat), rownames(terminal))
    expect_equal(unname(rowSums(mat)), rep(1, nrow(mat)), tolerance = 1e-12)
  }
  for (type in unique(h$subtype)) {
    expect_equal(fit$type.fraction[, type],
                 rowSums(terminal[, h$subtype == type, drop = FALSE]), tolerance = 1e-14)
  }
  for (lineage in unique(h$major_lineage)) {
    expect_equal(fit$lineage.fraction[, lineage],
                 rowSums(terminal[, h$major_lineage == lineage, drop = FALSE]), tolerance = 1e-14)
  }
}

test_that("standard outputs follow subtype, type, lineage in every hierarchy mode", {
  x <- output_fixture()
  local_mocked_bindings(load_immucellai2_reference = function() x$ref,
                        load_immucellai2_markers = function() rownames(x$ref))
  for (mode in c("tcell", "flat")) {
    for (unknown in c(FALSE, TRUE)) {
      fit <- run_immucellai2(x$bulk, hierarchy.mode = mode,
                            add.unknown = unknown, n.iter = 3)
      expect_fraction_levels(fit)
      expect_equal(ncol(fit$subtype.fraction), 5L + as.integer(unknown))
      for (field in c("subtype.fraction", "type.fraction", "lineage.fraction")) {
        expect_identical("UNKNOWN" %in% colnames(fit[[field]]), unknown)
      }
      if (mode == "tcell" && !unknown) {
        expect_equal(ncol(fit$type.fraction), 4L)
        expect_equal(ncol(fit$lineage.fraction), 3L)
      }
    }
  }
})

test_that("all matrix inference and auxiliary estimators use the same output names", {
  x <- output_fixture()
  for (method in c("vb", "mcmc")) {
    fit <- deconvolve_bulk_matrix(x$bulk, x$ref, hierarchy.mode = "tcell",
                                 inference.method = method, n.iter = 4, burn.in = 2)
    expect_fraction_levels(fit)
  }
  expect_fraction_levels(deconvolve_bulk_fast(x$bulk, x$ref, hierarchy.mode = "tcell",
                                             marker.n = 10, n.iter = 3))
  expect_fraction_levels(deconvolve_marker_score(x$bulk, x$ref, marker.n = 10))
})

test_that("exported fraction filenames round-trip at the correct hierarchy level", {
  x <- output_fixture()
  local_mocked_bindings(load_immucellai2_reference = function() x$ref,
                        load_immucellai2_markers = function() rownames(x$ref))
  fit <- run_immucellai2(x$bulk, n.iter = 3)
  out <- tempfile("output-names-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_true(write_deconvolution_outputs(fit, out))
  expect_setequal(list.files(out), c("subtype_fraction.txt", "type_fraction.txt",
                                     "lineage_fraction.txt", "hierarchy_used.txt"))
  for (level in c("subtype", "type", "lineage")) {
    actual <- as.matrix(read.table(file.path(out, paste0(level, "_fraction.txt")),
                                   header = TRUE, sep = "\t", row.names = 1,
                                   check.names = FALSE))
    expect_equal(actual, fit[[paste0(level, ".fraction")]], tolerance = 1e-14)
  }
})

test_that("legacy fits and directories cannot silently mix hierarchy levels", {
  x <- output_fixture()
  local_mocked_bindings(load_immucellai2_reference = function() x$ref,
                        load_immucellai2_markers = function() rownames(x$ref))
  fit <- run_immucellai2(x$bulk, n.iter = 3)
  old <- fit
  current <- c("subtype.fraction", "type.fraction", "lineage.fraction")
  legacy <- c("state.fraction", "subtype.fraction", "major.fraction")
  names(old)[match(current, names(old))] <- legacy
  out <- tempfile("legacy-output-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_error(write_deconvolution_outputs(old, out), "Legacy result names")
  expect_false(dir.exists(out))
  names(old)[match(legacy, names(old))] <- current
  expect_identical(old, fit)
  expect_true(write_deconvolution_outputs(old, out))
  for (filename in c("state_fraction.txt", "major_lineage_fraction.txt")) {
    legacy.path <- file.path(out, filename)
    writeLines("old output", legacy.path)
    before <- tools::md5sum(list.files(out, full.names = TRUE))
    expect_error(write_deconvolution_outputs(fit, out), "Use a new output directory")
    expect_identical(tools::md5sum(names(before)), before)
    unlink(legacy.path)
  }
})
