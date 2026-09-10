# Paths below are resolved relative to the repository root.
config <- list(
  data_dir = "benchmarks/fig2/data",
  output_dir = "benchmarks/output/fig2",
  # Optional named paths, e.g. c(GSE146771 = "path/to/simulation.txt").
  bulk_files = NULL,
  # NULL uses the atlas and marker list packaged with ImmuCellAI 2.0.
  reference_file = NULL,
  immune_gene_file = NULL,
  library_paths = character(0),
  n_cores = 8L,
  sample_limit = NA_integer_,
  chain_length = 600L,
  burn_in = 500L,
  cibersort_signature_genes = 1000L,
  run_immucellai2 = TRUE,
  run_bayesprism = TRUE,
  run_dwls = TRUE,
  run_music = TRUE,
  run_citmic = TRUE,
  run_cibersort = TRUE,
  run_immucellai = TRUE
)
