legacy_root <- getOption("Immu2.legacy_root", "C:/个人文件总库/02工作学习/R")
archive <- file.path(legacy_root, "Fig2/Newman_Monaco_7tools_validation/RNA_Sieve_newman_pbmcs")
ids <- c(ImmuCellAI2 = "ImmuneHier_flat_VB",
  BayesPrism = "BayesPrism_first_state_chain600_burn500", CIBERSORT = "CIBERSORT_default_LM22",
  CITMIC = "CITMIC_native", DWLS = "DWLS_native_ref_top1000", ImmuCellAI = "ImmuCellAI_native",
  MuSiC_basic = "MuSiC_native_music.basic")
relative <- setNames(paste0(ids, "_fraction.txt"), names(ids))
list(kind = "external", dataset = "RNA_Sieve_newman_pbmcs", archive = archive,
  bulk = file.path(archive, "standardized/expression.csv"),
  truth = file.path(archive, "standardized/true_proportions.csv"),
  reference = file.path(legacy_root, "reference_53celltypesTPM20260518.txt"),
  markers = file.path(legacy_root, "MarkerUsedDeconvolution.txt"),
  immunehier_mode = "flat", immunehier_add_unknown = FALSE,
  immunehier_method_name = "ImmuneHier_flat_VB",
  predictions = setNames(file.path(archive, "seven_tools_compare", relative), names(relative)),
  prediction_relative = relative, engine = "B/seven_tools.R",
  target_alias = c("Neutrophil" = "Neutrophils", "Monocyte" = "Monocytes",
    "CD8 T cell" = "CD8T", "CD4 T cell" = "CD4T", "B cell" = "B", "NK cell" = "NK"),
  plot_data = file.path(project, "Data/Fig3/B/pie_plot_data.csv"),
  note = "Historical flat configuration; UNKNOWN FALSE; native LM22 uses QN=TRUE here.")
