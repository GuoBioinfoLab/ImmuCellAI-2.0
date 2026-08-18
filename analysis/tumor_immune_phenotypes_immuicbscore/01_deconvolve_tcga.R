args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/tumor_immune_phenotypes_immuicbscore/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "ImmuCellAI2.0"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
tcga_expression_file <- assert_file(tcga_expression_file, "TCGA TPM matrix")
tcga_clinical_file <- assert_file(tcga_clinical_file, "TCGA clinical table")

clinical <- read_table_auto(tcga_clinical_file)
required <- c("case_submitter_id", "project_id")
if (!all(required %in% names(clinical))) {
  stop("TCGA clinical table must contain: ", paste(required, collapse = ", "))
}
case_project <- unique(clinical[clinical$project_id %in% tcga_projects, required])
case_to_project <- setNames(case_project$project_id, case_project$case_submitter_id)

header <- names(data.table::fread(tcga_expression_file, nrows = 0L, check.names = FALSE))
if (length(header) < 2L) stop("TCGA TPM matrix must contain genes in column 1 and samples thereafter.")
sample_ids <- header[-1L]
case_ids <- substr(sample_ids, 1L, 12L)
sample_type <- substr(sample_ids, 14L, 15L)
keep <- case_ids %in% names(case_to_project) & sample_type %in% c("01", "02", "03", "05", "06")
selected <- sample_ids[keep]
if (!length(selected)) stop("No TCGA samples matched the configured projects and tumor sample types.")

sample_table <- data.frame(
  sample = selected,
  case_submitter_id = substr(selected, 1L, 12L),
  project_id = unname(case_to_project[substr(selected, 1L, 12L)]),
  sample_type_code = substr(selected, 14L, 15L),
  stringsAsFactors = FALSE
)
data.table::fwrite(sample_table, file.path(output_dir, "tcga_selected_samples.tsv"), sep = "\t")

selected_file <- file.path(output_dir, "TCGA_selected_TPM.tsv")
dt <- data.table::fread(
  tcga_expression_file,
  select = c(header[1L], selected),
  check.names = FALSE,
  showProgress = TRUE
)
data.table::fwrite(dt, selected_file, sep = "\t", quote = FALSE)
rm(dt)
gc()

fit <- run_standard_deconvolution(
  expression_file = selected_file,
  output_dir = file.path(output_dir, "tcga_deconvolution"),
  n_cores = n_cores
)
message("TCGA deconvolution complete: ", nrow(fit$state.fraction), " samples.")
