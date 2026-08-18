repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- file.path(repo_root, "analysis", "tumor_immune_phenotypes_immuicbscore", "data")
output_dir <- file.path(repo_root, "analysis", "tumor_immune_phenotypes_immuicbscore", "results")

tcga_expression_file <- file.path(data_dir, "TCGA_TPM_symbol.txt")
tcga_clinical_file <- file.path(data_dir, "TCGA_clinical.tsv")
tcga_survival_file <- file.path(data_dir, "TCGA_survival.tsv")
tcga_fraction_file <- file.path(output_dir, "tcga_deconvolution", "state_fraction.txt")
tcga_cluster_file <- file.path(output_dir, "consensus_clustering", "consensus_class_k4.tsv")

icb_fraction_file <- file.path(data_dir, "icb_immune_fractions.tsv")
icb_clinical_file <- file.path(data_dir, "icb_clinical.tsv")
icb_cv_output_dir <- file.path(output_dir, "immuicbscore_fivefold_cv")

n_cores <- 8L
consensus_reps <- 1000L
consensus_seed <- 123L

tcga_projects <- c(
  "TCGA-SKCM", "TCGA-LUAD", "TCGA-LUSC", "TCGA-KIRC", "TCGA-BLCA",
  "TCGA-HNSC", "TCGA-STAD", "TCGA-ESCA", "TCGA-LIHC", "TCGA-CESC",
  "TCGA-UCEC", "TCGA-COAD", "TCGA-READ", "TCGA-MESO"
)

selected_icb_studies <- c(
  "anti-PD1_SRP070710",
  "anti-PD1_SRP230414",
  "anti-PD1_SRP351936",
  "anti-PD1_ERP105482",
  "anti-PD1-anti-CTLA4_ERP105482",
  "anti-PD1_ERP107734",
  "anti-PD1_ERP117672",
  "anti-CTLA4-to-anti-PD1_SRP417444"
)
