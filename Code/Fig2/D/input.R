cohort <- file.path(project, "GSE162025")
archive <- file.path(cohort, "02_seven_tools/selected_CD4_127_Mono_125_11types")
simulation <- file.path(cohort, "01_simulation/generated_selected_CD4_127_Mono_125_11types")
relative <- setNames(file.path("predictions", paste0(pipeline_methods, ".rds")), pipeline_methods)
list(kind = "npc", dataset = "GSE162025", archive = archive, cohort = cohort, simulation = simulation,
  bulk = file.path(simulation, "GSE162025_simulated_expression_440.tsv.gz"),
  pairs = file.path(archive, "metrics/sample_level_truth_prediction_primary.tsv"),
  mapping = file.path(archive, "input/locked_mapping_11types.tsv"),
  locked_protocol = file.path(archive, "input/locked_protocol.json"),
  predictions = setNames(file.path(archive, relative), names(relative)), prediction_relative = relative,
  engine = "D/seven_tools.R", plot_data = file.path(project, "Data/Fig2/D/pie_plot_data.tsv"),
  note = "Fixed CD4T127/Monocyte125 blocks. Evaluate all-cell truth; no immune-conditional renormalization.")
