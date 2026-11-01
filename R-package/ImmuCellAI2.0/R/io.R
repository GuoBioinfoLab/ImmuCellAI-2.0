#' Read an expression matrix
#'
#' The first column is interpreted as the gene identifier. Duplicate gene
#' identifiers are collapsed by their arithmetic mean. Values must be numeric;
#' missing and non-finite values are replaced by zero with a warning, whereas
#' negative values are rejected because the deconvolution model expects
#' non-negative abundances.
#'
#' @param file Tab-delimited expression matrix with a header row.
#' @return A numeric genes-by-samples matrix.
#' @export
read_expression_matrix <- function(file) {
  if (length(file) != 1L || is.na(file) || !nzchar(file) || !file.exists(file)) {
    stop("File not found: ", file)
  }

  x <- utils::read.table(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    row.names = NULL
  )
  if (ncol(x) < 2L) {
    stop("Expression file must contain one gene column and at least one sample column.")
  }

  genes <- trimws(as.character(x[[1L]]))
  if (anyNA(genes) || any(!nzchar(genes))) {
    stop("The first column contains missing or empty gene identifiers.")
  }

  values <- x[-1L]
  numeric.values <- lapply(values, function(v) suppressWarnings(as.numeric(v)))
  bad.numeric <- vapply(seq_along(values), function(i) {
    any(is.na(numeric.values[[i]]) & !is.na(values[[i]]) & nzchar(trimws(as.character(values[[i]]))))
  }, logical(1))
  if (any(bad.numeric)) {
    stop(
      "Non-numeric expression values were found in sample column(s): ",
      paste(names(values)[bad.numeric], collapse = ", ")
    )
  }

  mat <- as.matrix(as.data.frame(numeric.values, check.names = FALSE))
  storage.mode(mat) <- "double"
  colnames(mat) <- names(values)

  nonfinite <- !is.finite(mat)
  if (any(nonfinite)) {
    warning(sum(nonfinite), " missing or non-finite expression value(s) were replaced by zero.")
    mat[nonfinite] <- 0
  }
  if (any(mat < 0)) {
    stop("Negative expression values are not supported. Supply a non-negative TPM-like matrix.")
  }

  if (anyDuplicated(genes)) {
    warning("Duplicate gene identifiers were collapsed by arithmetic mean.")
    mat <- rowsum(mat, group = genes, reorder = FALSE) /
      as.numeric(table(factor(genes, levels = unique(genes))))
    rownames(mat) <- unique(genes)
  } else {
    rownames(mat) <- genes
  }

  if (anyDuplicated(colnames(mat))) {
    stop("Sample names must be unique.")
  }
  mat
}

#' Write deconvolution outputs
#'
#' Writes subtype_fraction.txt (terminal populations), type_fraction.txt
#' (middle-level groups), lineage_fraction.txt (top-level branches), and
#' hierarchy_used.txt. Legacy result objects and directories containing legacy
#' fraction filenames are rejected to prevent mixing hierarchy levels.
#'
#' @param result Result from deconvolve_bulk_matrix().
#' @param out.dir Output directory.
#' @export
write_deconvolution_outputs <- function(result, out.dir) {
  if (any(c("state.fraction", "major.fraction") %in% names(result))) {
    stop("Legacy result names detected. Migrate the result as documented in NEWS for v0.2.0 before exporting.")
  }
  legacy.files <- file.path(out.dir, c("state_fraction.txt", "major_lineage_fraction.txt"))
  if (any(file.exists(legacy.files))) {
    stop("Legacy output files detected. Use a new output directory to avoid mixing hierarchy levels.")
  }
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
  if (!is.null(result$subtype.fraction)) {
    utils::write.table(result$subtype.fraction, file.path(out.dir, "subtype_fraction.txt"),
                       sep = "\t", quote = FALSE, col.names = NA)
  }
  if (!is.null(result$type.fraction)) {
    utils::write.table(result$type.fraction, file.path(out.dir, "type_fraction.txt"),
                       sep = "\t", quote = FALSE, col.names = NA)
  }
  if (!is.null(result$lineage.fraction)) {
    utils::write.table(result$lineage.fraction, file.path(out.dir, "lineage_fraction.txt"),
                       sep = "\t", quote = FALSE, col.names = NA)
  }
  if (!is.null(result$hierarchy)) {
    utils::write.table(result$hierarchy, file.path(out.dir, "hierarchy_used.txt"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(TRUE)
}
