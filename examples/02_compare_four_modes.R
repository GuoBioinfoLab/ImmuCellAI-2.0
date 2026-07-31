suppressPackageStartupMessages(library(ImmuneHierDeconv))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript examples/02_compare_four_modes.R bulk_TPM.txt [output_dir] [n_cores]")
}

bulk <- read_expression_matrix(args[1])
output_dir <- if (length(args) >= 2L) args[2] else "ImmuCellAI2_four_modes"
n_cores <- if (length(args) >= 3L) as.integer(args[3]) else 1L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

design <- data.frame(
  mode = c("tcell_only", "tcell_only", "flat", "flat"),
  unknown = c(FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(design))) {
  label <- paste0(design$mode[i], if (design$unknown[i]) "_UNKNOWN" else "")
  message("Running ", label)
  fit <- run_immucellai2(
    bulk = bulk,
    hierarchy.mode = design$mode[i],
    add.unknown = design$unknown[i],
    n.cores = n_cores,
    seed = 123
  )
  write_deconvolution_outputs(fit, file.path(output_dir, label))
  saveRDS(fit$settings, file.path(output_dir, label, "settings.rds"))
}

utils::write.table(
  design,
  file.path(output_dir, "mode_design.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
