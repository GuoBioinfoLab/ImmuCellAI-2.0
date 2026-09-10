# Assemble the current Figure 3 from the original vector ggplot objects.
source(file.path(Sys.getenv("IMMUCELLAI2_FIG3_HOME", "benchmarks/fig3"), "paths.R"), local = TRUE)

if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package patchwork is required.")
library(patchwork)

panel_env <- new.env(parent = globalenv())
sys.source(
  fig3_script("plots/make_gse107011_newman_monaco_runtime_combined.R"),
  envir = panel_env
)
summary_env <- new.env(parent = globalenv())
sys.source(
  fig3_script("summary/plot_6datasets_7tools_mean_pearson.R"),
  envir = summary_env
)

tag_theme <- ggplot2::theme(
  plot.tag = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
  plot.tag.position = c(0, 1)
)

p_a <- panel_env$p_gse +
  ggplot2::facet_grid(
    . ~ Dataset,
    labeller = ggplot2::as_labeller(c(GSE107011 = "Xu's Cohort (PBMC; GSE107019)"))
  ) +
  ggplot2::labs(tag = "A") + tag_theme
p_b <- panel_env$p_newman +
  ggplot2::facet_grid(
    . ~ Dataset,
    labeller = ggplot2::as_labeller(c("Newman PBMC" = "Alizadeh's Cohort (PBMC)"))
  ) +
  ggplot2::labs(tag = "B") + tag_theme
p_c <- panel_env$p_monaco +
  ggplot2::facet_grid(
    . ~ Dataset,
    labeller = ggplot2::as_labeller(c("Monaco PBMC" = "Larbi's Cohort (PBMC)"))
  ) +
  ggplot2::labs(tag = "C") + tag_theme
p_d <- panel_env$p_runtime + ggplot2::labs(tag = "D") + tag_theme

summary_labels <- c(
  GSE146771 = "CRC (GSE146771)",
  GSE164522 = "CRC (GSE164522)",
  GSE176078 = "BRCA (GSE176078)",
  GSE107011 = "PBMC (Xu; GSE107019)",
  "Newman PBMC" = "PBMC (Alizadeh)",
  "Monaco PBMC" = "PBMC (Larbi)"
)
p_e <- summary_env$p +
  ggplot2::facet_wrap(
    ~ Dataset, ncol = 6,
    labeller = ggplot2::as_labeller(summary_labels)
  ) +
  ggplot2::labs(tag = "E") + tag_theme

middle_row <- p_b + p_c + p_d +
  patchwork::plot_layout(ncol = 3, widths = c(1.05, 1, 0.92))
figure3 <- p_a / middle_row / p_e +
  patchwork::plot_layout(heights = c(1.12, 1, 0.82))

out_dir <- file.path(fig3_work_dir(), "Figure3_complete")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(out_dir, "Figure3_PBMC_validation_benchmark_complete.pdf")
out_png <- file.path(out_dir, "Figure3_PBMC_validation_benchmark_complete.png")
ggplot2::ggsave(out_pdf, figure3, width = 11.69, height = 11.55,
                units = "in", bg = "white", limitsize = FALSE)
ggplot2::ggsave(out_png, figure3, width = 11.69, height = 11.55,
                units = "in", dpi = 500, bg = "white", limitsize = FALSE)
message("Saved complete Figure 3 PDF: ", out_pdf)
message("Saved complete Figure 3 PNG: ", out_png)
