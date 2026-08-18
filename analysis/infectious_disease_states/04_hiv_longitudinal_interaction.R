args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/infectious_disease_states/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2", "readxl"))

fraction <- read_fraction_matrix(fraction_file)
clinical <- read_table_auto(hiv_longitudinal_file)
required <- c("Run", "Patient_ID", "grouping", "outcome")
if (!all(required %in% names(clinical))) stop("Longitudinal metadata requires: ", paste(required, collapse = ", "))

grouping_to_time <- function(x) dplyr::case_when(
  grepl("_0$", x) ~ "D0", grepl("_1$", x) ~ "D1",
  grepl("_4i?$", x) ~ "D4", grepl("_8$", x) ~ "D8", TRUE ~ NA_character_
)
normalize_outcome <- function(x) dplyr::case_when(
  tolower(trimws(x)) %in% c("died", "dead") ~ "Died",
  tolower(trimws(x)) %in% c("survived", "survive") ~ "Survived", TRUE ~ as.character(x)
)
time_day <- c(D0 = 0, D1 = 1, D4 = 4, D8 = 8)
main_cells <- c("Tc", "cNK", "Neutrophil", "MDSC")
supplement_cells <- c("Th1", "CD4Tn", "CD8Tn", "Tr1")
cells <- c(main_cells, supplement_cells)
missing <- setdiff(cells, colnames(fraction))
if (length(missing)) stop("Missing longitudinal cells: ", paste(missing, collapse = ", "))

dat <- data.frame(Run = rownames(fraction), fraction[, cells, drop = FALSE], check.names = FALSE) |>
  dplyr::left_join(clinical, by = "Run") |>
  dplyr::mutate(
    status = factor(normalize_outcome(outcome), levels = c("Died", "Survived")),
    Time = factor(grouping_to_time(grouping), levels = names(time_day)),
    TimeDay = unname(time_day[as.character(Time)])
  ) |>
  dplyr::filter(!is.na(status), !is.na(Time)) |>
  tidyr::pivot_longer(dplyr::all_of(cells), names_to = "CellType", values_to = "Fraction")

fit_one <- function(x) {
  repeated <- any(table(x$Patient_ID) > 1L)
  if (repeated && requireNamespace("lmerTest", quietly = TRUE)) {
    fit <- tryCatch(
      lmerTest::lmer(Fraction ~ status * TimeDay + (1 | Patient_ID), data = x, REML = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      tab <- as.data.frame(summary(fit)$coefficients)
      term <- "statusSurvived:TimeDay"
      if (term %in% rownames(tab)) {
        return(data.frame(model = "mixed-effects", estimate = tab[term, "Estimate"], se = tab[term, "Std. Error"], p = tab[term, "Pr(>|t|)"]))
      }
    }
  }
  fit <- lm(Fraction ~ status * TimeDay, data = x)
  tab <- as.data.frame(summary(fit)$coefficients)
  term <- "statusSurvived:TimeDay"
  if (!term %in% rownames(tab)) {
    return(data.frame(model = "not-estimable", estimate = NA_real_, se = NA_real_, p = NA_real_))
  }
  data.frame(model = "linear", estimate = tab[term, "Estimate"], se = tab[term, "Std. Error"], p = tab[term, "Pr(>|t|)"])
}
stats <- dat |>
  dplyr::group_by(CellType) |>
  dplyr::group_modify(~fit_one(.x)) |>
  dplyr::ungroup() |>
  dplyr::mutate(FDR = p.adjust(p, "BH"))
summary <- dat |>
  dplyr::group_by(CellType, status, Time, TimeDay) |>
  dplyr::summarise(mean = mean(Fraction), sem = sd(Fraction) / sqrt(dplyr::n()), .groups = "drop")
data.table::fwrite(dat, file.path(output_dir, "HIV_longitudinal_sample_data.tsv"), sep = "\t")
data.table::fwrite(stats, file.path(output_dir, "HIV_status_by_time_interaction.tsv"), sep = "\t")

make_plot <- function(selected, stem) {
  x <- summary |> dplyr::filter(CellType %in% selected) |>
    dplyr::mutate(CellType = factor(CellType, levels = selected))
  p <- ggplot2::ggplot(x, ggplot2::aes(Time, mean, color = status, group = status)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sem, ymax = mean + sem), width = 0.12) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(~CellType, nrow = 1, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(Died = "#A3C8DC", Survived = "#EA5D2D")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(linetype = "dashed", color = "grey80"), aspect.ratio = 1) +
    ggplot2::labs(x = "Time", y = "Mean cell fraction", color = NULL)
  save_plot_pair(p, file.path(output_dir, stem), 12, 3.2)
}
make_plot(main_cells, "HIV_longitudinal_main_cells")
make_plot(supplement_cells, "HIV_longitudinal_supplement_cells")
