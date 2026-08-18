# Age-associated immune-cell remodeling

This module reproduces the healthy peripheral-blood age case study.

1. `01_prepare_and_deconvolve.R` merges the two TPM matrices, aligns age
   metadata, and runs ImmuCellAI 2.0.
2. `02_plot_age_cell_abundance.R` creates the manuscript abundance panels.
3. `03_plot_age_pattern_heatmap.R` summarizes age-group means, scales each cell
   state across age groups, and draws the ordered heatmap.
4. `04_plot_major_lineage_composition.R` draws the stacked major-lineage plot.

Copy `config.example.R` to `config.R`, edit paths, and run:

```bash
Rscript analysis/age_associated_remodeling/run_all.R \
  analysis/age_associated_remodeling/config.R
```

Display trimming in panel plots does not remove samples from statistical input
or exported tables. The unmodified inferred fractions remain in the
deconvolution output.
