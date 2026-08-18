args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/infectious_disease_states/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2", "ggpubr", "gghalves"))

fraction <- read_fraction_matrix(fraction_file)
meta <- read_table_auto(tb_metadata_file)
if (!all(c("Run", "tb_status") %in% names(meta))) stop("TB metadata requires Run and tb_status.")
cells <- c("cDC1", "Tc", "cMo", "MDSC")
missing <- setdiff(cells, colnames(fraction))
if (length(missing)) stop("Missing TB plot cells: ", paste(missing, collapse = ", "))
dat <- data.frame(Run = rownames(fraction), fraction[, cells, drop = FALSE], check.names = FALSE) |>
  dplyr::left_join(meta, by = "Run") |>
  dplyr::filter(tb_status %in% c("Active TB", "Latent TB")) |>
  dplyr::mutate(tb_status = factor(tb_status, levels = c("Active TB", "Latent TB"))) |>
  tidyr::pivot_longer(dplyr::all_of(cells), names_to = "CellType", values_to = "Fraction") |>
  dplyr::mutate(CellType = factor(CellType, levels = cells))
data.table::fwrite(dat, file.path(output_dir, "TB_active_vs_latent_plot_data.tsv"), sep = "\t")

p <- ggplot2::ggplot(dat, ggplot2::aes(tb_status, Fraction, fill = tb_status)) +
  gghalves::geom_half_violin(side = "r", trim = TRUE, alpha = 0.9, color = "black") +
  ggplot2::geom_boxplot(width = 0.28, outlier.shape = NA, color = "black", alpha = 0.9) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", size = 5) +
  ggplot2::facet_wrap(~CellType, nrow = 1, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c("Active TB" = "#D9534F", "Latent TB" = "#4E79A7")) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_line(linetype = "dashed", color = "grey80"), aspect.ratio = 1) +
  ggplot2::labs(x = NULL, y = "Cell fraction")
save_plot_pair(p, file.path(output_dir, "TB_Active_vs_Latent_selected_cells"), 12, 3.2)

stats <- dat |>
  dplyr::group_by(CellType) |>
  dplyr::summarise(
    n_active = sum(tb_status == "Active TB"), n_latent = sum(tb_status == "Latent TB"),
    mean_active = mean(Fraction[tb_status == "Active TB"]),
    mean_latent = mean(Fraction[tb_status == "Latent TB"]),
    p_value = wilcox.test(Fraction ~ tb_status)$p.value, .groups = "drop"
  ) |>
  dplyr::mutate(FDR = p.adjust(p_value, "BH"))
data.table::fwrite(stats, file.path(output_dir, "TB_active_vs_latent_statistics.tsv"), sep = "\t")
