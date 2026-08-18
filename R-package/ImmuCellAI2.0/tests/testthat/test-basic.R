test_that("basic deconvolution runs", {
  set.seed(1)
  ref <- matrix(rexp(100 * 4), nrow = 100, ncol = 4)
  rownames(ref) <- paste0("Gene", 1:100)
  colnames(ref) <- paste0("Cell", 1:4)
  theta <- c(0.7, 0.1, 0.1, 0.1)
  bulk <- ref %*% theta
  bulk <- matrix(bulk, ncol = 1)
  rownames(bulk) <- rownames(ref)
  colnames(bulk) <- "S1"
  hierarchy <- data.frame(
    state_name = colnames(ref),
    subtype = colnames(ref),
    major_lineage = "Immune",
    stringsAsFactors = FALSE
  )
  res <- deconvolve_bulk_matrix(bulk, ref, hierarchy, n.iter = 5, burn.in = 2, n.cores = 1)
  expect_equal(nrow(res$state.fraction), 1)
  expect_equal(ncol(res$state.fraction), 4)
})

test_that("variational Bayesian deconvolution runs", {
  set.seed(3)
  ref <- matrix(rexp(120 * 4), nrow = 120, ncol = 4)
  rownames(ref) <- paste0("Gene", 1:120)
  colnames(ref) <- paste0("Cell", 1:4)
  theta <- c(0.55, 0.25, 0.15, 0.05)
  bulk <- matrix(ref %*% theta, ncol = 1)
  rownames(bulk) <- rownames(ref)
  colnames(bulk) <- "S1"
  hierarchy <- data.frame(
    state_name = colnames(ref),
    subtype = c("A", "A", "B", "B"),
    major_lineage = c("Immune1", "Immune1", "Immune2", "Immune2"),
    stringsAsFactors = FALSE
  )
  res <- deconvolve_bulk_matrix(
    bulk,
    ref,
    hierarchy,
    inference.method = "vb",
    n.iter = 20,
    n.cores = 1
  )
  expect_equal(nrow(res$state.fraction), 1)
  expect_equal(ncol(res$state.fraction), 4)
  expect_equal(round(sum(res$state.fraction[1, ]), 8), 1)
  expect_equal(res$inference.method, "vb")
})

test_that("flat variational Bayesian deconvolution runs", {
  set.seed(5)
  ref <- matrix(rexp(120 * 4), nrow = 120, ncol = 4)
  rownames(ref) <- paste0("Gene", 1:120)
  colnames(ref) <- paste0("Cell", 1:4)
  theta <- c(0.55, 0.25, 0.15, 0.05)
  bulk <- matrix(ref %*% theta, ncol = 1)
  rownames(bulk) <- rownames(ref)
  colnames(bulk) <- "S1"
  hierarchy <- data.frame(
    state_name = colnames(ref),
    subtype = c("A", "A", "B", "B"),
    major_lineage = c("Immune1", "Immune1", "Immune2", "Immune2"),
    stringsAsFactors = FALSE
  )
  res <- deconvolve_bulk_matrix(
    bulk,
    ref,
    hierarchy,
    inference.method = "flat_vb",
    n.iter = 20,
    n.cores = 1
  )
  expect_equal(nrow(res$state.fraction), 1)
  expect_equal(ncol(res$state.fraction), 4)
  expect_equal(round(sum(res$state.fraction[1, ]), 8), 1)
  expect_equal(res$inference.method, "flat_vb")
})

test_that("tcell hierarchy leaves non-T states independent", {
  h <- create_tcell_53_hierarchy(c("CD4Tcm", "CD4Tem", "CD8Tcm", "BGC", "cMo"))
  expect_equal(h$major_lineage[h$state_name == "CD4Tcm"], "CD4T")
  expect_equal(h$subtype[h$state_name == "CD4Tem"], "CD4 Memory T")
  expect_equal(h$major_lineage[h$state_name == "CD8Tcm"], "CD8T")
  expect_equal(h$major_lineage[h$state_name == "BGC"], "BGC")
  expect_equal(h$subtype[h$state_name == "cMo"], "cMo")
})


test_that("tcell_only remains a backward-compatible alias", {
  states <- c("CD4Tcm", "CD4Tem", "CD8Tcm", "BGC", "cMo")
  canonical <- create_default_53_hierarchy(states, hierarchy.mode = "tcell")
  legacy <- create_default_53_hierarchy(states, hierarchy.mode = "tcell_only")
  expect_equal(legacy, canonical)
  expect_equal(immucellai2_defaults()$hierarchy.mode, "tcell")
})
test_that("fast marker score and alias evaluation run", {
  set.seed(2)
  ref <- matrix(rexp(120 * 5), nrow = 120, ncol = 5)
  rownames(ref) <- paste0("Gene", 1:120)
  colnames(ref) <- c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm", "cMo")
  theta <- c(0.1, 0.2, 0.3, 0.1, 0.3)
  bulk <- matrix(ref %*% theta, ncol = 1)
  rownames(bulk) <- rownames(ref)
  colnames(bulk) <- "CD4Tmemory_1"

  hierarchy <- create_default_53_hierarchy(colnames(ref), hierarchy.mode = "hybrid")
  res <- deconvolve_marker_score(bulk, ref, hierarchy = hierarchy, marker.n = 10)
  expect_equal(nrow(res$state.fraction), 1)
  expect_equal(ncol(res$state.fraction), 5)
  expect_equal(round(sum(res$state.fraction[1, ]), 8), 1)

  sample.info <- extract_sample_truth_from_names(colnames(bulk))
  ev <- evaluate_prediction_matrix(
    res$state.fraction,
    sample.info,
    target.mapper = map_simulation_target_to_states
  )
  expect_true("CD4Tcm+CD4Tem+CD4Temra+CD4Trm" %in% ev$sample.info$EvalTarget)
  expect_false(is.na(ev$sample.info$Predicted[1]))
})

test_that("expression reader preserves names and collapses duplicate genes", {
  input <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "gene\tS1\tS2",
      "CD3D\t1\t3",
      "MS4A1\t2\t4",
      "CD3D\t5\t7"
    ),
    input
  )

  expect_warning(x <- read_expression_matrix(input), "Duplicate gene")
  expect_equal(rownames(x), c("CD3D", "MS4A1"))
  expect_equal(colnames(x), c("S1", "S2"))
  expect_equal(unname(x["CD3D", ]), c(3, 5))
})

test_that("expression reader rejects negative and non-numeric values", {
  negative <- tempfile(fileext = ".txt")
  writeLines(c("gene\tS1", "CD3D\t-1"), negative)
  expect_error(read_expression_matrix(negative), "Negative expression")

  nonnumeric <- tempfile(fileext = ".txt")
  writeLines(c("gene\tS1", "CD3D\tnot_a_number"), nonnumeric)
  expect_error(read_expression_matrix(nonnumeric), "Non-numeric expression")
})


test_that("expression reader handles an unlabeled gene column", {
  input <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "CellA\tCellB",
      "Gene1\t1\t2",
      "Gene2\t3\t4"
    ),
    input
  )

  x <- read_expression_matrix(input)
  expect_equal(rownames(x), c("Gene1", "Gene2"))
  expect_equal(colnames(x), c("CellA", "CellB"))
  expect_equal(dim(x), c(2L, 2L))
})
