# Tumor immune phenotypes and ImmuICBscore

This module contains configurable analyses and 47 full local source scripts for
the TCGA immune-phenotype and immunotherapy case study. See the
[source preparation guide](../REPRODUCING_FIGURES.md) before running archived files.

## Figure 4 panel map

The paths below are relative to `source_archive/` unless stated otherwise.

| Panel | Analysis and plotting implementations |
|---|---|
| A | `run_tcga_immunotherapy_responsive_immucellai2_marker5000.R` prepares the selected TCGA fractions; `ConsensusClusterPlus.r` and `Fig4AB_Code.R` contain clustering variants; `ComplexHeatmap.r` draws the ordered four-IP heatmap. `additional/plot_tcga_k4_immune_heatmap.R` provides a separate K=4 heatmap export. |
| B | `analyze_tcga_cluster_score_survival_association.R` contains the KM and Cox analyses; the top-level `03_plot_phenotypes_and_survival.R` is the configurable four-IP survival entry point. |
| C-D | `plot_icb334_pre_on_response_cell_fraction.R`: cDC1, Tc, CD8Tem and MBC, pre/on-treatment within response groups and R/NR comparisons on treatment. |
| E | `fig4_01_unified_icb_score_method.R`: CLR, scaling, class-balanced probability forest, prediction; `fig4_04_plot_unified_icb_score_method_diagram.R` and its `icb_development/` variant draw method diagrams from the bundled background assets. Final publication layout is separately composed. |
| F-G | `analyze_unified_icb_score_in_tcga_clusters.R` trains the descriptive TCGA scoring model; `comprehensive_cluster4_unified_icb_score_association.R` writes score/cluster summaries and quartiles; `fig4_02_plot_unified_icb_score_cluster_and_roc.R` plots the cluster-score distributions and stacked quartiles. |
| H | `icb_development/fig4_07_fivefold_cv_unified_icb_score.R` through `fig4_10_search_best_8th_cohort_cv_roc.R` preserve the historical CV and cohort-selection workflow. The last script exports `best_8cohort_pooled_5fold_predictions.txt`, summary, ROC coordinates, and ROC plots for the 334-sample result. Use the top-level `04_immuicbscore_pooled_fivefold_cv.R` for the already fixed eight-cohort evaluation. |

The ROC drawn by the earlier `fig4_02` script uses LOSO outputs and is **not**
the final 334-sample panel H. Likewise, `Fig4AB_Code.R` retains a K=5 development
section and interactive objects; do not mistake that section for the final
four-phenotype assignment. Use saved K=4 labels or the configurable K=4 workflow.

## Historical execution and inputs

TCGA analysis requires gene-by-sample TPM, project metadata, saved K=4 sample
labels, and survival metadata (historical `survival_pan_all.RData` or the
documented tabular equivalent). ICB analysis requires sample-by-53 fractions and
matching response, cohort and treatment-time metadata. The historical external
ICB deconvolution/preparation is in `validate_C3_C4_score_external_immunotherapy.R`.
Downstream inputs are listed at the beginning of each script; preserve this
dependency order: fractions and labels, fitted score/summary tables, then plots.

For F-G, run `analyze_unified_icb_score_in_tcga_clusters.R`, then
`comprehensive_cluster4_unified_icb_score_association.R`, then `fig4_02...`.
The last script also requires the historical LOSO files for its older ROC;
the configurable `05_associate_score_with_tcga_phenotypes.R` avoids that coupling.
TCGA scores are descriptive model outputs, not observed treatment responses.

The archived C-D script uses two-sided, unpaired Wilcoxon tests, even for
pre/on-treatment groups. Its P values should not be described as paired tests.
The core forest uses 1,200 trees, `mtry = floor(sqrt(p))`, minimum node size 8,
probability output and inverse-frequency class weights. The reference `p=53`
gives `mtry=7`; exploratory scripts may drop constant features or vary settings.

The remaining `optimize_*`, signature, cluster-trend, stage and marker-comparison
scripts are retained for transparency, not as additional required main panels.
Source hashes and dependencies are in [the manifest](../source_manifest.json).

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
