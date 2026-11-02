cohort <- file.path(project, "real_bulk_flow_validation_20261011/additional_cohorts_round2/BG_blood")
archive <- file.path(cohort, "seven_tools_50iter_selected26")
relative <- setNames(file.path("predictions", paste0(pipeline_methods, ".rds")), pipeline_methods)
list(kind = "bg", dataset = "BostonGene blood", archive = archive, cohort = cohort,
  bulk = file.path(cohort, "bulk_linear.tsv.gz"), truth = file.path(cohort, "truth_published.tsv"),
  pairs = file.path(archive, "metrics/sample_level_truth_prediction.tsv"),
  mapping = file.path(archive, "input/locked_mapping.tsv"),
  locked_protocol = file.path(archive, "input/locked_protocol.json"),
  baseline = file.path(cohort, "results/ImmuCellAI2_53states.tsv"),
  predictions = setNames(file.path(archive, relative), names(relative)), prediction_relative = relative,
  engine = "A/seven_tools.R", plot_data = file.path(project, "Data/Fig3/A/pie_plot_data.csv"),
  note = "45 RNA samples; 26 displayed targets. The original runner reuses the independently verified 50-iteration ImmuCellAI2 baseline.")
