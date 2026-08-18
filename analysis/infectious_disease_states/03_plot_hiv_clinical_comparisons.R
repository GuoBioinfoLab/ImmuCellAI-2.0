args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/infectious_disease_states/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2", "ggpubr", "gghalves", "patchwork"))

fraction <- read_fraction_matrix(fraction_file)
meta <- read_table_auto(hiv_metadata_file)
if (!all(c("Run", "outcome", "treatment") %in% names(meta))) {
  stop("HIV metadata requires Run, outcome, and treatment.")
}
normalize_outcome <- function(x) dplyr::case_when(
  tolower(trimws(x)) %in% c("died", "dead") ~ "Died",
  tolower(trimws(x)) %in% c("survived", "survive") ~ "Survived",
  tolower(trimws(x)) %in% c("died-iris", "dead-iris") ~ "Died-IRIS",
  tolower(trimws(x)) %in% c("survived-iris", "survive-iris") ~ "Survived-IRIS",
  TRUE ~ as.character(x)
)
meta$outcome <- normalize_outcome(meta$outcome)
base <- data.frame(Run = rownames(fraction), fraction, check.names = FALSE) |>
  dplyr::left_join(meta, by = "Run")

plot_comparison <- function(groups, cells, label, colors) {
  missing <- setdiff(cells, colnames(fraction))
  if (length(missing)) stop("Missing cells for ", label, ": ", paste(missing, collapse = ", "))
  dat <- base |>
    dplyr::filter(outcome %in% groups) |>
    dplyr::mutate(Group = factor(outcome, levels = groups)) |>
    tidyr::pivot_longer(dplyr::all_of(cells), names_to = "CellType", values_to = "Fraction") |>
    dplyr::mutate(CellType = factor(CellType, levels = cells))
  p <- ggplot2::ggplot(dat, ggplot2::aes(Group, Fraction, fill = Group)) +
    gghalves::geom_half_violin(side = "r", trim = TRUE, color = "black") +
    ggplot2::geom_boxplot(width = 0.28, outlier.shape = NA, color = "black") +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif") +
    ggplot2::facet_wrap(~CellType, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = colors) + ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_line(linetype = "dashed", color = "grey80"), aspect.ratio = 1) +
    ggplot2::labs(x = NULL, y = "Cell fraction")
  save_plot_pair(p, file.path(output_dir, label), 12, 3.2)
  data.table::fwrite(dat, file.path(output_dir, paste0(label, "_data.tsv")), sep = "\t")
}

outcome_cells <- c("MDSC", "Neutrophil", "Th1", "Tc")
plot_comparison(
  c("Died", "Survived"), outcome_cells, "HIV_Died_vs_Survived",
  c(Died = "#E69F3A", Survived = "#59A14F")
)
plot_comparison(
  c("Died-IRIS", "Survived-IRIS"), outcome_cells, "HIV_DiedIRIS_vs_SurvivedIRIS",
  c(`Died-IRIS` = "#E69F3A", `Survived-IRIS` = "#59A14F")
)

art_cells <- c("cNK", "CD8Temra", "ncMo", "MDSC")
art <- base |>
  dplyr::filter(treatment %in% c("Early ART", "Deferred ART")) |>
  dplyr::mutate(Group = factor(treatment, levels = c("Early ART", "Deferred ART"))) |>
  tidyr::pivot_longer(dplyr::all_of(art_cells), names_to = "CellType", values_to = "Fraction") |>
  dplyr::mutate(CellType = factor(CellType, levels = art_cells))
p_art <- ggplot2::ggplot(art, ggplot2::aes(Group, Fraction, fill = Group)) +
  gghalves::geom_half_violin(side = "r", trim = TRUE, color = "black") +
  ggplot2::geom_boxplot(width = 0.28, outlier.shape = NA, color = "black") +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif") +
  ggplot2::facet_wrap(~CellType, nrow = 1, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(`Early ART` = "#59A14F", `Deferred ART` = "#EDC948")) +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
  ggplot2::labs(x = NULL, y = "Cell fraction")
save_plot_pair(p_art, file.path(output_dir, "HIV_Early_vs_Deferred_ART"), 12, 3.2)
data.table::fwrite(art, file.path(output_dir, "HIV_Early_vs_Deferred_ART_data.tsv"), sep = "\t")
