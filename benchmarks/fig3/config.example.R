config <- list(
  work_dir = "benchmarks/output/fig3",
  # NULL uses work_dir/inputs for the Xu expression, truth and ID-mapping files.
  input_dir = NULL,
  fig2_results_dir = "benchmarks/output/fig2",
  reference_file = NULL,
  immune_gene_file = NULL,
  library_paths = character(0),
  n_cores = 8L,
  # TRUE explicitly permits the original runtime plotter's archived constants.
  use_archived_timing = FALSE
)
