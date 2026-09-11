message("Running Fig4 Unified ICB Score integration workflow...")

code_dir <- "<LOCAL_R_ROOT>/Fig4/Unified_ICB_score_cluster_ICB_integration/code"

scripts <- c(
  "fig4_01_unified_icb_score_method.R",
  "fig4_02_plot_unified_icb_score_cluster_and_roc.R",
  "fig4_04_plot_unified_icb_score_method_diagram.R",
  "fig4_03_extended_tcga_cluster_score_analysis.R",
  "fig4_05_fixed_validation_unified_icb_score.R",
  "fig4_06_optimize_train_validation_auc.R"
)

for (script in scripts) {
  script_path <- file.path(code_dir, script)
  if (!file.exists(script_path)) {
    stop("Missing script: ", script_path)
  }
  message("Sourcing: ", script_path)
  source(script_path, encoding = "UTF-8")
}

message("Fig4 workflow complete.")
