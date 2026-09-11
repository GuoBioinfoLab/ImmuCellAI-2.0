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

frac_file <- file.path(out_dir, "ImmuCellAI2_PRJNA683803_state_fraction_sample_by_celltype.txt")
clinical_file <- file.path(fig6_dir, "PRJNA683803.xlsx")

main_cells <- c("Tc", "cNK", "Neutrophil", "MDSC")

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
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fit_one <- function(dat) {
  dat <- dat %>%
    mutate(
      status = factor(status, levels = c("Died", "Survived")),
      TimeDay = as.numeric(TimeDay)
    )
  fit <- lm(Abundance ~ status * TimeDay, data = dat)
  coef_tab <- as.data.frame(summary(fit)$coefficients)
  term <- "statusSurvived:TimeDay"
  data.frame(
    Model = "linear_model_lm",
    Term = term,
    Estimate = coef_tab[term, "Estimate"],
    StdError = coef_tab[term, "Std. Error"],
    Pvalue = coef_tab[term, "Pr(>|t|)"],
    N = nrow(dat),
    stringsAsFactors = FALSE
  )
}

make_interaction_plot <- function(plot_cells, long_dat, stats_dat, suffix, title_text,
                                  width = 12.5, height = 6.0) {
  plot_summary <- long_dat %>%
    filter(CellType %in% plot_cells) %>%
    mutate(CellType = factor(CellType, levels = plot_cells)) %>%
    group_by(CellType, status, Time, TimeDay) %>%
    summarise(
      N = n(),
      Mean = mean(Abundance, na.rm = TRUE),
      SEM = sd(Abundance, na.rm = TRUE) / sqrt(N),
      .groups = "drop"
    )

  plot_points <- long_dat %>%
    filter(CellType %in% plot_cells) %>%
    mutate(CellType = factor(CellType, levels = plot_cells))

  labels <- stats_dat %>%
    filter(CellType %in% plot_cells) %>%
    mutate(CellType = factor(CellType, levels = plot_cells)) %>%
    transmute(
      CellType,
      label = paste0("Pint(D0-D1) = ", format_p(Pvalue), "\nFDR = ", format_p(FDR))
    )

  y_pos <- plot_points %>%
    group_by(CellType) %>%
    summarise(y = max(Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = y + 0.08 * abs(y))

  labels <- labels %>% left_join(y_pos, by = "CellType")

  phy_cols <- c("Died" = "#A3C8DC", "Survived" = "#EA5D2D")
  p <- ggplot(plot_summary, aes(x = Time, y = Mean, group = status, color = status)) +
    geom_point(
      data = plot_points,
      aes(x = Time, y = Abundance, color = status),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.08, height = 0),
      alpha = 0.18,
      size = 1.0,
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
      size = 3.0,
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
         p, width = width, height = height, bg = "white")
  ggsave(file.path(out_dir, paste0("Figure6_PRJNA683803_status_time_interaction_", suffix, ".png")),
         p, width = width, height = height, dpi = 600, bg = "white")
  invisible(p)
}

message("Reading data...")
frac <- fread(frac_file, data.table = FALSE, check.names = FALSE)
colnames(frac)[1] <- "Run"
clinical <- readxl::read_excel(clinical_file) %>% as.data.frame()
all_cells <- setdiff(colnames(frac), "Run")

dat <- frac %>%
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
  select(Run, Patient_ID, status, Time, TimeDay, all_of(all_cells)) %>%
  pivot_longer(cols = all_of(all_cells), names_to = "CellType", values_to = "Abundance") %>%
  mutate(Abundance = as.numeric(Abundance))

early_long <- long %>% filter(Time %in% c("D0", "D1"))

early_stats <- early_long %>%
  group_by(CellType) %>%
  group_modify(~ fit_one(.x)) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(Pvalue, method = "BH")) %>%
  arrange(FDR, Pvalue)

sig8_cells <- early_stats %>%
  filter(!CellType %in% main_cells) %>%
  arrange(FDR, Pvalue) %>%
  slice_head(n = 8) %>%
  pull(CellType)

write.table(early_stats, file.path(out_dir, "PRJNA683803_status_time_interaction_all53_early_D0D1_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(CellType = sig8_cells),
            file.path(out_dir, "PRJNA683803_status_time_interaction_supplement_sig8_cells.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

make_interaction_plot(
  sig8_cells,
  long,
  early_stats,
  "supplement_sig8_cells",
  "Died vs Survived early status-by-time interaction"
)

cat("Selected supplement cells:\n")
print(early_stats %>% filter(CellType %in% sig8_cells) %>%
        select(CellType, Estimate, StdError, Pvalue, FDR, N))
