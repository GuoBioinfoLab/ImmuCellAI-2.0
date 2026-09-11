options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

fig6_dir <- "<LOCAL_R_ROOT>/Fig6"
out_dir <- file.path(fig6_dir, "ImmuCellAI2_PRJNA683803_AIDS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

state_meta_file <- file.path(out_dir, "ImmuCellAI2_PRJNA683803_state_fraction_with_metadata.txt")
metadata_file <- file.path(fig6_dir, "PRJNA683803_combined.csv")

read_optional_grouping <- function(fig6_dir) {
  candidates <- file.path(
    fig6_dir,
    c("PRJNA683803_grouping.csv", "PRJNA683803_grouping.txt",
      "PRJNA683803_grouping.tsv", "PRJNA683803_grouping.xlsx",
      "PRJNA683803.xlsx",
      "PRJNA683803_combined_with_grouping.csv",
      "PRJNA683803_combined_with_grouping.txt",
      "PRJNA683803_combined_with_grouping.xlsx")
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) return(NULL)
  f <- candidates[1]
  if (grepl("\\.xlsx$", f, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Found grouping xlsx but package 'readxl' is not installed: ", f)
    }
    dat <- readxl::read_excel(f) %>% as.data.frame()
  } else {
    dat <- fread(f, data.table = FALSE, check.names = FALSE)
  }
  dat
}

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

message("Reading ImmuCellAI2 PRJNA683803 state fractions with metadata...")
state_df <- fread(state_meta_file, data.table = FALSE, check.names = FALSE)

metadata <- fread(metadata_file, data.table = FALSE, check.names = FALSE)
metadata$outcome <- normalise_outcome(metadata$outcome)

if (!"grouping" %in% colnames(metadata)) {
  grouping_df <- read_optional_grouping(fig6_dir)
  if (!is.null(grouping_df)) {
    run_col <- intersect(c("Run", "run", "Sample", "sample", "SampleID", "sample_id"), colnames(grouping_df))[1]
    grouping_col <- intersect(c("grouping", "Grouping", "group", "Group"), colnames(grouping_df))[1]
    if (is.na(run_col) || is.na(grouping_col)) {
      stop("The grouping file must contain a sample column such as 'Run' and a grouping column such as 'grouping'.")
    }
    metadata <- metadata %>%
      left_join(grouping_df %>% select(Run = all_of(run_col), grouping = all_of(grouping_col)), by = "Run")
  }
}

if (!"grouping" %in% colnames(metadata) || all(is.na(metadata$grouping))) {
  stop(
    "Cannot draw the p4 time-course plot because PRJNA683803_combined.csv has no 'grouping' column.\n",
    "The original p4 code requires sample-level groupings such as ED-I_0, DS-I_1, ES-N_4, DS-N_8 to map samples to D0/D1/D4/D8.\n",
    "Please add a 'grouping' column to PRJNA683803_combined.csv, or place a file named PRJNA683803_grouping.csv/txt/xlsx in Fig6 with columns Run and grouping."
  )
}

metadata <- metadata %>%
  mutate(
    outcome = normalise_outcome(outcome),
    Time = grouping_to_time(grouping),
    Time = factor(Time, levels = c("D0", "D1", "D4", "D8")),
    status = factor(outcome, levels = c("Died", "Survived"))
  )

merged <- state_df %>%
  select(-any_of(c("BioProject", "Disease", "outcome", "treatment", "grouping", "Time", "status"))) %>%
  left_join(metadata %>% select(Run, outcome, status, grouping, Time), by = "Run") %>%
  filter(status %in% c("Died", "Survived"), !is.na(Time))

main_cells <- c("Tc", "cNK", "Neutrophil", "MDSC")
supplement_cells <- c("CD4Tn", "Th1", "Tr1", "CD8Tn")
selected_cells <- unique(c(main_cells, supplement_cells))
missing_cells <- setdiff(selected_cells, colnames(merged))
if (length(missing_cells) > 0L) {
  stop("Missing selected ImmuCellAI2 cells: ", paste(missing_cells, collapse = ", "))
}

phy_cols <- c("Died" = "#A3C8DC", "Survived" = "#EA5D2D")

make_timecourse_plot <- function(cells, table_prefix, figure_prefix) {
  plot_long <- merged %>%
    select(Run, status, Time, all_of(cells)) %>%
    pivot_longer(cols = all_of(cells), names_to = "CellType", values_to = "Abundance") %>%
    mutate(
      Abundance = as.numeric(Abundance),
      CellType = factor(CellType, levels = cells)
    )

  summary_long <- plot_long %>%
    group_by(Time, status, CellType) %>%
    summarise(
      Mean = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    )

  write.table(plot_long, file.path(out_dir, paste0(table_prefix, "_sample_level_long.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary_long, file.path(out_dir, paste0(table_prefix, "_group_mean_long.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  p <- ggplot(summary_long, aes(x = Time, y = Mean, group = status, color = status)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    facet_wrap(~ CellType, ncol = 4, scales = "free_y") +
    labs(x = "Time", y = "Ratios", title = "Cell Type Ratios by Status") +
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
    ) +
    scale_colour_manual(values = phy_cols, drop = FALSE)

  ggsave(file.path(out_dir, paste0(figure_prefix, ".pdf")),
         p, width = 12.5, height = 3.2, bg = "white")
  ggsave(file.path(out_dir, paste0(figure_prefix, ".png")),
         p, width = 12.5, height = 3.2, dpi = 600, bg = "white")

  invisible(list(plot = p, sample_level = plot_long, group_mean = summary_long))
}

main_result <- make_timecourse_plot(
  cells = main_cells,
  table_prefix = "PRJNA683803_timecourse_p4_main",
  figure_prefix = "Figure6_PRJNA683803_timecourse_p4_Died_vs_Survived"
)

supplement_result <- make_timecourse_plot(
  cells = supplement_cells,
  table_prefix = "PRJNA683803_timecourse_p4_supplement",
  figure_prefix = "Figure6_PRJNA683803_timecourse_p4_supplement_Died_vs_Survived"
)

cat("Done. Time-course p4 plot written to:\n", out_dir, "\n", sep = "")
