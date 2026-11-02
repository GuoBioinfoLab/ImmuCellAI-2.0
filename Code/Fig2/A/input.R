# Frozen 40-point target gradients; no simulation or seed selection here.
archive <- file.path(project, "Fig2_GSE164522")
ids <- c(ImmuCellAI2 = "ImmuneHier_tcell_only_VB",
  BayesPrism = "BayesPrism_first_state_chain600_burn500", CIBERSORT = "CIBERSORT_LM22",
  CITMIC = "CITMIC_native", DWLS = "DWLS_weighted_lm", ImmuCellAI = "ImmuCellAI_native",
  MuSiC_basic = "MuSiC_basic")
relative <- setNames(paste0(ids, "_fraction.txt"), names(ids))
relative["CIBERSORT"] <- "corrections/CIBERSORT_LM22/CIBERSORT_LM22_predictions.rds"
predictions <- setNames(file.path(archive, "historical/predictions", paste0(ids, "_fraction.txt")), names(ids))
predictions["CIBERSORT"] <- file.path(archive, relative["CIBERSORT"])
list(kind = "legacy", dataset = "GSE164522", archive = archive,
  bulk = file.path(archive, "input/bulk_original.txt"),
  truth = file.path(archive, "input/target_truth_manifest.tsv"),
  reference = file.path(archive, "references/reference_53celltypesTPM20260518.txt"),
  markers = file.path(archive, "references/MarkerUsedDeconvolution.txt"),
  predictions = predictions, prediction_relative = relative,
  engine = "A/seven_tools.R",
  plot_data = file.path(project, "Data/Fig2/A/pie_plot_data.tsv"),
  note = "Six historical sample-level prediction files are missing; only corrected LM22 predictions are available.")
