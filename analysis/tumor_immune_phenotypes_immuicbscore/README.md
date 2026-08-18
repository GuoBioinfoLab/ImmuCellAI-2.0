# Tumor immune phenotypes and ImmuICBscore

This module reproduces the manuscript case study on TCGA immune phenotypes and
immunotherapy response assessment.

## Workflow

1. `01_deconvolve_tcga.R`: select configured TCGA tumor projects and estimate
   53 immune-state fractions with ImmuCellAI 2.0.
2. `02_consensus_clustering.R`: consensus PAM clustering and K=4 immune
   phenotype assignment.
3. `03_plot_phenotypes_and_survival.R`: phenotype heatmap and Kaplan-Meier
   overall-survival analysis.
4. `04_immuicbscore_pooled_fivefold_cv.R`: pooled stratified five-fold
   cross-validation of the CLR-transformed, class-balanced ranger model.
5. `05_associate_score_with_tcga_phenotypes.R`: fit the final descriptive model
   on the selected ICB cohort set, score TCGA samples, and plot phenotype-score
   associations.
6. `06_plot_icb_response_cell_fractions.R`: pre/on-treatment and responder/
   non-responder cell-fraction comparisons.

Copy `config.example.R` to `config.R`, edit paths, then run each numbered script
or:

```bash
Rscript analysis/tumor_immune_phenotypes_immuicbscore/run_all.R \
  analysis/tumor_immune_phenotypes_immuicbscore/config.R
```

Required packages include `ImmuCellAI2.0`, `data.table`,
`ConsensusClusterPlus`, `ComplexHeatmap`, `circlize`, `survival`, `survminer`,
`ranger`, `dplyr`, `tidyr`, `ggplot2`, `ggpubr`, and `scales`.

The pooled AUC is a cross-validated estimate after exploratory cohort selection,
not a prospective treatment-homogeneous validation result.
