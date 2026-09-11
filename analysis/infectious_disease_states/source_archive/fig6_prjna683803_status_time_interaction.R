options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
})

fig6_dir <- "<LOCAL_R_ROOT>/Fig6"
out_dir <- file.path(fig6_dir, "ImmuCellAI2_PRJNA683803_AIDS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

frac_file <- file.path(out_dir, "ImmuCellAI2_PRJNA683803_state_fraction_sample_by_celltype.txt")
clinical_file <- file.path(fig6_dir, "PRJNA683803.xlsx")

main_cells <- c("Tc", "cNK", "Neutrophil", "MDSC")
supplement_cells <- c("Tc", "cNK", "Neutrophil", "MDSC", "Th1", "CD4Tn", "CD8Tn", "Tr1")
selected_cells <- unique(c(main_cells, supplement_cells))

normalise_outcome <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    tolower(x) %in% c("died", "dead") ~ "Died",
    tolower(x) %in% c("survived", "survive") ~ "Survived",
    tolower(x) %in% c("died-iris", "dead-iris") ~ "Died-IRIS",
    tolower(x) %in% c("survived-iris", "survive-iris") ~ "Survived-IRIS",
    TRUE ~ x
  )
}

grouping_to_time <- function(grouping) {
  grouping <- trimws(as.character(grouping))
  dplyr::case_when(
    grouping %in% c("ED-I_0", "DS-I_0", "ES-I_0", "ES-N_0",
                    "DD-N_0", "ED-N_0", "DS-N_0") ~ "D0",
    grouping %in% c("ED-I_1", "DS-I_1", "ES-I_1", "ES-N_1",
                    "DD-N_1", "ED-N_1", "DS-N_1") ~ "D1",
    grouping %in% c("ED-I_4i", "DS-I_4i", "ES-I_4i", "ES-N_4", "DS-N_4") ~ "D4",
    grouping %in% c("DS-I_8", "ES-I_8", "ES-N_8", "DS-N_8") ~ "D8",
    TRUE ~ NA_character_
  )
}

time_to_day <- function(time) {
  dplyr::case_when(
    time == "D0" ~ 0,
    time == "D1" ~ 1,
    time == "D4" ~ 4,
    time == "D8" ~ 8,
    TRUE ~ NA_real_
  )
}

format_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fit_interaction_model <- function(dat) {
  dat <- dat %>%
    mutate(
      status = factor(status, levels = c("Died", "Survived")),
      TimeDay = as.numeric(TimeDay),
      Patient_ID = factor(Patient_ID)
    )

  n_patient <- n_distinct(dat$Patient_ID)
  n_repeated <- dat %>%
    count(Patient_ID) %>%
    filter(n > 1) %>%
    nrow()

  if (requireNamespace("lmerTest", quietly = TRUE) &&
      requireNamespace("lme4", quietly = TRUE) &&
      n_patient >= 3 && n_repeated >= 2) {
    fit <- tryCatch(
      lmerTest::lmer(Abundance ~ status * TimeDay + (1 | Patient_ID), data = dat, REML = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      coef_tab <- as.data.frame(summary(fit)$coefficients)
      term <- "statusSurvived:TimeDay"
      if (term %in% rownames(coef_tab)) {
        return(data.frame(
          Model = "mixed_effects_lmer",
          Term = term,
          Estimate = coef_tab[term, "Estimate"],
          StdError = coef_tab[term, "Std. Error"],
          Pvalue = coef_tab[term, "Pr(>|t|)"],
          N = nrow(dat),
          N_patient = n_patient,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  fit <- lm(Abundance ~ status * TimeDay, data = dat)
  coef_tab <- as.data.frame(summary(fit)$coefficients)
  term <- "statusSurvived:TimeDay"
  data.frame(
    Model = "linear_model_lm",
    Term = term,
    Estimate = if (term %in% rownames(coef_tab)) coef_tab[term, "Estimate"] else NA_real_,
    StdError = if (term %in% rownames(coef_tab)) coef_tab[term, "Std. Error"] else NA_real_,
    Pvalue = if (term %in% rownames(coef_tab)) coef_tab[term, "Pr(>|t|)"] else NA_real_,
    N = nrow(dat),
    N_patient = n_patient,
    stringsAsFactors = FALSE
  )
}

message("Reading ImmuCellAI2 fractions and clinical grouping...")
frac <- fread(frac_file, data.table = FALSE, check.names = FALSE)
colnames(frac)[1] <- "Run"
clinical <- readxl::read_excel(clinical_file) %>% as.data.frame()

if (!all(c("Run", "Patient_ID", "grouping", "outcome") %in% colnames(clinical))) {
  stop("PRJNA683803.xlsx must contain Run, Patient_ID, grouping and outcome columns.")
}

missing_cells <- setdiff(selected_cells, colnames(frac))
if (length(missing_cells) > 0L) {
  stop("Missing selected ImmuCellAI2 cells: ", paste(missing_cells, collapse = ", "))
}

dat <- frac %>%
  select(Run, all_of(selected_cells)) %>%
  left_join(
    clinical %>%
      transmute(
        Run = as.character(Run),
        Patient_ID = as.character(Patient_ID),
        outcome = normalise_outcome(outcome),
        grouping = as.character(grouping),
        Time = grouping_to_time(grouping)
      ),
    by = "Run"
  ) %>%
  filter(outcome %in% c("Died", "Survived"), !is.na(Time)) %>%
  mutate(
    status = factor(outcome, levels = c("Died", "Survived")),
    Time = factor(Time, levels = c("D0", "D1", "D4", "D8")),
    TimeDay = time_to_day(as.character(Time))
  )

long <- dat %>%
  select(Run, Patient_ID, status, Time, TimeDay, all_of(selected_cells)) %>%
  pivot_longer(cols = all_of(selected_cells), names_to = "CellType", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(Abundance),
    CellType = factor(CellType, levels = selected_cells)
  )

early_long <- long %>% filter(Time %in% c("D0", "D1"))
full_long <- long

early_stats <- early_long %>%
  group_by(CellType) %>%
  group_modify(~ fit_interaction_model(.x)) %>%
  ungroup() %>%
  mutate(
    Analysis = "early_D0_D1_common_time_window",
    FDR = p.adjust(Pvalue, method = "BH")
  )

full_stats <- full_long %>%
  group_by(CellType) %>%
  group_modify(~ fit_interaction_model(.x)) %>%
  ungroup() %>%
  mutate(
    Analysis = "full_D0_D8_unbalanced_descriptive",
    FDR = p.adjust(Pvalue, method = "BH")
  )

interaction_stats <- bind_rows(early_stats, full_stats) %>%
  select(Analysis, CellType, Model, Term, Estimate, StdError, Pvalue, FDR, N, N_patient)

summary_by_time <- full_long %>%
  group_by(CellType, status, Time, TimeDay) %>%
  summarise(
    N = n(),
    Mean = mean(Abundance, na.rm = TRUE),
    SEM = sd(Abundance, na.rm = TRUE) / sqrt(N),
    .groups = "drop"
  )

write.table(long, file.path(out_dir, "PRJNA683803_status_time_interaction_sample_level_long.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_by_time, file.path(out_dir, "PRJNA683803_status_time_interaction_group_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(interaction_stats, file.path(out_dir, "PRJNA683803_status_time_interaction_model_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

phy_cols <- c("Died" = "#A3C8DC", "Survived" = "#EA5D2D")

make_interaction_plot <- function(cells, suffix, title_text, height = 3.2) {
  plot_summary <- summary_by_time %>%
    filter(CellType %in% cells) %>%
    mutate(CellType = factor(CellType, levels = cells))

  plot_points <- full_long %>%
    filter(CellType %in% cells) %>%
    mutate(CellType = factor(CellType, levels = cells))

  labels <- early_stats %>%
    filter(CellType %in% cells) %>%
    transmute(
      CellType = factor(CellType, levels = cells),
      label = paste0("Pint(D0-D1) = ", format_p(Pvalue), "\nFDR = ", format_p(FDR))
    )

  y_pos <- plot_points %>%
    group_by(CellType) %>%
    summarise(y = max(Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = y + 0.08 * abs(y))

  labels <- labels %>% left_join(y_pos, by = "CellType")

  p <- ggplot(plot_summary, aes(x = Time, y = Mean, group = status, color = status)) +
    geom_point(
      data = plot_points,
      aes(x = Time, y = Abundance, color = status),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.08, height = 0),
      alpha = 0.18,
      size = 1.1,
      show.legend = FALSE
    ) +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                  width = 0.12, linewidth = 0.55,
                  position = position_dodge(width = 0.08)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.3) +
    geom_text(
      data = labels,
      aes(x = "D0", y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.1,
      color = "black"
    ) +
    facet_wrap(~ CellType, ncol = 4, scales = "free_y") +
    labs(x = "Time", y = "Cell fraction", title = title_text) +
    scale_colour_manual(values = phy_cols, drop = FALSE) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.75, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.05),
      aspect.ratio = 1,
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 14, color = "black", hjust = 0),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )

  ggsave(file.path(out_dir, paste0("Figure6_PRJNA683803_status_time_interaction_", suffix, ".pdf")),
         p, width = 12.5, height = height, bg = "white")
  ggsave(file.path(out_dir, paste0("Figure6_PRJNA683803_status_time_interaction_", suffix, ".png")),
         p, width = 12.5, height = height, dpi = 600, bg = "white")
  invisible(p)
}

make_effect_plot <- function(stats_df) {
  plot_df <- stats_df %>%
    filter(Analysis == "early_D0_D1_common_time_window") %>%
    mutate(
      CellType = factor(CellType, levels = rev(selected_cells)),
      Lower = Estimate - 1.96 * StdError,
      Upper = Estimate + 1.96 * StdError,
      Panel = ifelse(as.character(CellType) %in% main_cells, "Main cells", "Supplement cells")
    )

  p <- ggplot(plot_df, aes(x = Estimate, y = CellType, color = Panel)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.18, linewidth = 0.65) +
    geom_point(size = 2.4) +
    labs(
      x = "Interaction beta: Survived x time day, D0-D1",
      y = NULL,
      title = "Status-by-time interaction effect"
    ) +
    scale_color_manual(values = c("Main cells" = "#EA5D2D", "Supplement cells" = "#4E79A7FF")) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.6, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.05),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 13, color = "black"),
      plot.title = element_text(size = 14, color = "black", hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(file.path(out_dir, "Figure6_PRJNA683803_status_time_interaction_effect_size.pdf"),
         p, width = 6.8, height = 4.6, bg = "white")
  ggsave(file.path(out_dir, "Figure6_PRJNA683803_status_time_interaction_effect_size.png"),
         p, width = 6.8, height = 4.6, dpi = 600, bg = "white")
  invisible(p)
}

make_interaction_plot(
  main_cells,
  "main_cells",
  "Died vs Survived status-by-time interaction"
)
make_interaction_plot(
  supplement_cells,
  "supplement_cells",
  "Died vs Survived early status-by-time interaction",
  height = 6.0
)
make_effect_plot(interaction_stats)

cat("Done. Interaction model outputs written to:\n", out_dir, "\n", sep = "")
