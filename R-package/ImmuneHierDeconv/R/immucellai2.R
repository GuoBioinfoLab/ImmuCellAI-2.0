#' Published ImmuCellAI 2.0 analysis settings
#'
#' Returns the settings used for the principal analyses in the ImmuCellAI 2.0
#' manuscript. The optional UNKNOWN residual is disabled in the published
#' workflow and can be enabled explicitly for sensitivity analysis.
#'
#' @return A named list of deconvolution settings.
#' @export
immucellai2_defaults <- function() {
  list(
    hierarchy.mode = "tcell_only",
    inference.method = "vb",
    add.unknown = FALSE,
    pseudo.depth = 1e5,
    n.iter = 50,
    vb.tol = 1e-6,
    alpha.major = 10,
    alpha.sub = 5,
    alpha.state = 1,
    seed = 123
  )
}

#' Load the ImmuCellAI 2.0 53-state reference atlas
#'
#' @param file Optional path to a tab-delimited reference matrix. The first
#'   column must contain gene symbols and the remaining columns must contain
#'   the 53 immune-state profiles.
#' @return A genes-by-cell-states numeric matrix.
#' @export
load_immucellai2_reference <- function(file = NULL) {
  if (is.null(file)) {
    file <- system.file(
      "extdata", "reference_53celltypesTPM20260518.txt",
      package = "ImmuneHierDeconv"
    )
  }
  if (!nzchar(file) || !file.exists(file)) {
    stop("ImmuCellAI 2.0 reference atlas was not found. Supply `file` explicitly.")
  }
  read_expression_matrix(file)
}

#' Load the ImmuCellAI 2.0 marker-gene panel
#'
#' @param file Optional path to a one-gene-per-row text file. A header named
#'   `gene` is accepted.
#' @return A character vector of unique gene symbols.
#' @export
load_immucellai2_markers <- function(file = NULL) {
  if (is.null(file)) {
    file <- system.file(
      "extdata", "MarkerUsedDeconvolution_5510.txt",
      package = "ImmuneHierDeconv"
    )
  }
  if (!nzchar(file) || !file.exists(file)) {
    stop("ImmuCellAI 2.0 marker-gene file was not found. Supply `file` explicitly.")
  }
  genes <- trimws(readLines(file, warn = FALSE))
  genes <- genes[nzchar(genes)]
  if (length(genes) > 0 && tolower(genes[1]) %in% c("gene", "genes", "symbol")) {
    genes <- genes[-1]
  }
  unique(genes)
}

#' Run the standard ImmuCellAI 2.0 workflow
#'
#' This is the recommended user-facing entry point. It intersects bulk RNA-seq
#' TPM values with the 53-state reference atlas and the published 5,510-gene
#' marker panel, then runs deterministic variational Bayesian inference with
#' targeted hierarchical refinement of CD4 and CD8 T-cell states.
#'
#' @param bulk A genes-by-samples numeric matrix or a path to a tab-delimited
#'   matrix whose first column contains gene symbols.
#' @param reference Optional genes-by-cell-states matrix or reference file.
#'   If omitted, the packaged 53-state atlas is used.
#' @param marker.genes Optional character vector or marker file. If omitted,
#'   the packaged 5,510-gene panel is used.
#' @param hierarchy.mode One of `tcell_only`, `flat`, `hierarchical`, or
#'   `hybrid`. The manuscript workflow uses `tcell_only`.
#' @param add.unknown Add an optional residual state for sensitivity analysis.
#'   The manuscript workflow uses `FALSE`.
#' @param n.cores Number of samples processed concurrently.
#' @param pseudo.depth Pseudo-count depth used after within-sample scaling.
#' @param n.iter Maximum number of deterministic VB updates.
#' @param vb.tol Convergence tolerance for the maximum fraction change.
#' @param alpha.major,alpha.sub,alpha.state Hierarchical Dirichlet prior
#'   concentrations.
#' @param seed Reproducibility seed.
#' @param verbose Print per-sample progress when running on one core.
#' @return The result from [deconvolve_bulk_matrix()] plus a `settings` field.
#' @export
run_immucellai2 <- function(
  bulk,
  reference = NULL,
  marker.genes = NULL,
  hierarchy.mode = "tcell_only",
  add.unknown = FALSE,
  n.cores = 1,
  pseudo.depth = 1e5,
  n.iter = 50,
  vb.tol = 1e-6,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  seed = 123,
  verbose = FALSE
) {
  if (is.character(bulk) && length(bulk) == 1L) {
    bulk <- read_expression_matrix(bulk)
  }
  if (is.null(reference)) {
    reference <- load_immucellai2_reference()
  } else if (is.character(reference) && length(reference) == 1L) {
    reference <- read_expression_matrix(reference)
  }
  if (is.null(marker.genes)) {
    marker.genes <- load_immucellai2_markers()
  } else if (is.character(marker.genes) && length(marker.genes) == 1L && file.exists(marker.genes)) {
    marker.genes <- load_immucellai2_markers(marker.genes)
  }

  bulk <- as.matrix(bulk)
  reference <- as.matrix(reference)
  if (is.null(rownames(bulk)) || is.null(colnames(bulk))) {
    stop("`bulk` must have gene row names and sample column names.")
  }
  if (is.null(rownames(reference)) || is.null(colnames(reference))) {
    stop("`reference` must have gene row names and cell-state column names.")
  }

  common <- intersect(marker.genes, intersect(rownames(reference), rownames(bulk)))
  if (length(common) < 1000L) {
    warning("Only ", length(common), " marker genes overlap bulk and reference matrices.")
  }
  if (length(common) == 0L) {
    stop("No marker genes overlap bulk and reference matrices. Check gene identifiers.")
  }

  bulk.use <- bulk[common, , drop = FALSE]
  reference.use <- reference[common, , drop = FALSE]
  hierarchy <- create_default_53_hierarchy(
    colnames(reference.use),
    hierarchy.mode = hierarchy.mode
  )

  result <- deconvolve_bulk_matrix(
    bulk.mat = bulk.use,
    reference = reference.use,
    hierarchy = hierarchy,
    hierarchy.mode = hierarchy.mode,
    inference.method = "vb",
    add.unknown = add.unknown,
    pseudo.depth = pseudo.depth,
    n.iter = n.iter,
    vb.tol = vb.tol,
    alpha.major = alpha.major,
    alpha.sub = alpha.sub,
    alpha.state = alpha.state,
    n.cores = n.cores,
    seed = seed,
    verbose = verbose
  )
  result$settings <- list(
    hierarchy.mode = hierarchy.mode,
    inference.method = "vb",
    add.unknown = add.unknown,
    pseudo.depth = pseudo.depth,
    n.iter = n.iter,
    vb.tol = vb.tol,
    alpha.major = alpha.major,
    alpha.sub = alpha.sub,
    alpha.state = alpha.state,
    n.cores = n.cores,
    seed = seed,
    marker.genes.available = length(marker.genes),
    marker.genes.used = length(common)
  )
  result
}
