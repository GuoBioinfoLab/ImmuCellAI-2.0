archive <- file.path(project, "GSE176078_14types_mean_profile_simulation",
                     "exploratory_Bmemory_cDC2_20261004")
relative <- setNames(file.path("predictions", paste0(pipeline_methods, "_full_predictions.rds")), pipeline_methods)
list(kind = "brca", dataset = "GSE176078", archive = archive,
  cohort = file.path(project, "Fig2_GSE176078"),
  selected_dir = file.path(archive, "selected_360"),
  bulk = file.path(archive, "selected_360/simulated_expression.tsv.gz"),
  pairs = file.path(archive, "metrics/sample_level_truth_prediction.tsv"),
  mapping = file.path(archive, "metrics/mapping_9targets_7tools.tsv"),
  design = file.path(archive, "locked_simulation_design.rds"),
  reference = file.path(project, "Fig2_GSE176078/references/reference_53celltypesTPM20260518.txt"),
  predictions = setNames(file.path(archive, relative), names(relative)), prediction_relative = relative,
  engine = "C/seven_tools.R", plot_data = file.path(project, "Data/Fig2/C/pie_plot_data.tsv"),
  note = "Frozen selected 360 samples. Full seven-tool replay only; no new selection or incremental cache mixing.")
