args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/tumor_immune_phenotypes_immuicbscore/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2", "ggpubr"))

fraction <- read_fraction_matrix(icb_fraction_file)
clinical <- read_table_auto(icb_clinical_file)
required <- c("Run", "ResponseBinary", "TreatmentTime")
if (!all(required %in% names(clinical))) {
  stop("For pre/on-treatment plots, clinical data require Run, ResponseBinary, and TreatmentTime columns.")
}
cells <- c("cDC1", "Tc", "CD8Tem", "MBC")
missing <- setdiff(cells, colnames(fraction))
if (length(missing)) stop("Missing cell states: ", paste(missing, collapse = ", "))
dat <- data.frame(Run = rownames(fraction), fraction[, cells, drop = FALSE], check.names = FALSE)
dat <- merge(dat, clinical, by = "Run")
dat$Response <- factor(ifelse(dat$ResponseBinary == 1, "R", "NR"), levels = c("NR", "R"))
dat$TreatmentTime <- factor(dat$TreatmentTime, levels = c("Pre", "On"))
long <- tidyr::pivot_longer(dat, dplyr::all_of(cells), names_to = "CellType", values_to = "Fraction")

p_on <- long |>
  dplyr::filter(TreatmentTime == "On") |>
  ggplot2::ggplot(ggplot2::aes(Response, Fraction, fill = Response)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.85) +
  ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white") +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif") +
  ggplot2::facet_wrap(~CellType, nrow = 1, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(NR = "#4E79A7", R = "#E15759")) +
  ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Cell fraction")
save_plot_pair(p_on, file.path(output_dir, "ICB_on_treatment_R_vs_NR_selected_cells"), 12, 3.2)

p_time <- ggplot2::ggplot(long, ggplot2::aes(TreatmentTime, Fraction, fill = TreatmentTime)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.85) +
  ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white") +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif") +
  ggplot2::facet_grid(Response ~ CellType, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(Pre = "#76B7B2", On = "#F28E2B")) +
  ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Cell fraction")
save_plot_pair(p_time, file.path(output_dir, "ICB_pre_on_within_response_groups"), 12, 6)
data.table::fwrite(long, file.path(output_dir, "ICB_selected_cell_fraction_plot_data.tsv"), sep = "\t")
