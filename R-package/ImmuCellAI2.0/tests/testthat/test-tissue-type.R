tissue_fixture <- function() {
  states <- c("CD4Tn", "CD4Tcm", "CD8Tn", "Bnaive", "TAM", "MDSC")
  ref <- matrix(1, 1201L, length(states),
    dimnames = list(paste0("Gene", seq_len(1201L)), states))
  for (j in seq_along(states)) ref[(200L * (j - 1L) + 1L):(200L * j), j] <- 15
  weights <- cbind(S1 = c(5, 4, 3, 2, 1, 6), S2 = c(2, 2, 1, 4, 6, 3))
  list(ref = ref, bulk = ref %*% weights, markers = rownames(ref)[seq_len(1200L)])
}

test_that("non-tumor removes reference profiles before standard inference", {
  d <- tissue_fixture()
  local_mocked_bindings(load_immucellai2_reference = function() d$ref,
                        load_immucellai2_markers = function() d$markers)

  expect_identical(immucellai2_defaults()$tissue.type, "tumor")
  default <- run_immucellai2(d$bulk)
  tumor <- run_immucellai2(d$bulk, tissue.type = "tumor")
  expect_identical(default, tumor)
  expect_setequal(colnames(tumor$subtype.fraction), colnames(d$ref))

  ref <- d$ref[d$markers, !colnames(d$ref) %in% c("TAM", "MDSC"), drop = FALSE]
  for (mode in c("tcell", "flat")) {
    fit <- run_immucellai2(d$bulk, hierarchy.mode = mode, tissue.type = "non-tumor")
    expected <- deconvolve_bulk_matrix(d$bulk[d$markers, ], ref,
      hierarchy.mode = mode, inference.method = "vb", n.iter = 50)
    for (key in c("subtype.fraction", "type.fraction", "lineage.fraction",
                  "hierarchy", "genes.used")) {
      expect_identical(fit[[key]], expected[[key]])
    }
    expect_setequal(colnames(fit$subtype.fraction), colnames(ref))
    expect_equal(unname(rowSums(fit$subtype.fraction)), c(1, 1), tolerance = 1e-12)
    expect_identical(fit$settings$tissue.type, "non-tumor")
    expect_identical(fit$settings$hierarchy.mode, mode)
    expect_identical(fit$settings$marker.genes.used, length(d$markers))
  }

  residual <- run_immucellai2(d$bulk, tissue.type = "non-tumor",
    add.unknown = TRUE, n.iter = 2L)
  expect_setequal(colnames(residual$subtype.fraction), c(colnames(ref), "UNKNOWN"))
  expect_equal(unname(rowSums(residual$subtype.fraction)), c(1, 1), tolerance = 1e-12)
})

test_that("tissue type is explicit and applies only to the bundled reference", {
  for (invalid in list(NULL, NA_character_, 1, "non", c("tumor", "non-tumor"))) {
    expect_error(run_immucellai2(NULL, tissue.type = invalid), "tissue.type")
  }
  d <- tissue_fixture()
  expect_error(run_immucellai2(d$bulk, reference = d$ref, tissue.type = "non-tumor"),
    "only to the bundled reference")
})
