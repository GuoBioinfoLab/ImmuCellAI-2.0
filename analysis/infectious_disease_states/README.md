# Immune-cell differences across infectious disease states

This module contains configurable workflows and seven full original scripts for
the tuberculosis and HIV-associated disease applications. See the
[source preparation guide](../REPRODUCING_FIGURES.md).

## Figure 6 panel map

| Panel | Original script under `source_archive/` |
|---|---|
| A | `fig6_active_latent_tb_immucellai2_plot.R`: expression/metadata alignment, deconvolution, active/latent TB comparisons and half-violin/boxplots for cDC1, Tc, cMo and MDSC. |
| B-D | `fig6_prjna683803_aids_immucellai2_analysis.R`: PRJNA683803 deconvolution, outcome/IRIS/ART comparisons, Wilcoxon tables and selected/all-cell plots. B and C use MDSC, Neutrophil, Th1 and Tc; D uses MDSC, CD4Tn, Th1 and gdT. |
| E | `fig6_prjna683803_aids_timecourse_p4.R`: maps grouping labels to D0/D1/D4/D8 and draws outcome-stratified means for Tc, cNK, Neutrophil and MDSC. |

Run the TB script independently for A. Run the PRJNA683803 main script before
the time-course script for B-E. Their output folders are
`ImmuCellAI2_active_latent_TB/` and `ImmuCellAI2_PRJNA683803_AIDS/` under
`<LOCAL_R_ROOT>/Fig6/`. Both use `exp.txt` but different metadata tables.
The time-course script additionally needs grouping labels from a supported
grouping table or `PRJNA683803.xlsx`; see [data/README.md](data/README.md).

The original PRJNA683803 group-comparison code applies within-group IQR filtering
before both the displayed distributions and Wilcoxon tests. Preserve its
unfiltered exports when evaluating sensitivity to this choice. Panel E is a
descriptive mean trajectory, not a subject-level matched longitudinal test.

Additional full implementations are retained:
`fig6_prjna683803_status_time_interaction.R` (status-by-time modeling),
`fig6_select_early_interaction_sig8.R` (exploratory early-time selection),
`fig6_check_tex_pb_tb.R` (additional TB checks), and `TuberAIDS_analysis.r`
(older metadata assembly and interactive analyses). The interaction scripts
are not the analysis generating the significance symbols in A-D.

## Configurable workflow

1. `01_prepare_and_deconvolve.R` aligns public run identifiers and runs
   ImmuCellAI 2.0 once for all included samples.
2. `02_plot_active_vs_latent_tb.R` compares cDC1, Tc, cMo, and MDSC between
   active and latent tuberculosis.
3. `03_plot_hiv_clinical_comparisons.R` compares survival, IRIS, and ART timing
   groups using the reported selected cell populations.
4. `04_hiv_longitudinal_interaction.R` plots D0/D1/D4/D8 trajectories and tests
   status-by-time interactions, using a patient random intercept when repeated
   observations support a mixed-effects model.

Copy `config.example.R` to `config.R`, edit paths, and run:

```bash
Rscript analysis/infectious_disease_states/run_all.R \
  analysis/infectious_disease_states/config.R
```

Required plotting/model packages include `data.table`, `dplyr`, `tidyr`,
`ggplot2`, `ggpubr`, `gghalves`, `patchwork`, `readxl`, `lme4`, and `lmerTest`.
