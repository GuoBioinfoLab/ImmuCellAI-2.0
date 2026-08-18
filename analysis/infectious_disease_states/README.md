# Immune-cell differences across infectious disease states

This module contains the complete deconvolution, statistical analysis, and
plotting workflow for the tuberculosis and HIV-associated disease applications.

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
