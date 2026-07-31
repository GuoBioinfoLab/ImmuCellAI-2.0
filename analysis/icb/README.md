# ImmuICBscore analysis

## Scope

ImmuICBscore is a probability score trained from the 53 relative immune
fractions returned by ImmuCellAI 2.0. It uses no gene-expression signatures,
IFN-gamma score, cytotoxicity score, checkpoint score, MHC score, or T-cell
inflamed GEP.

## Model

For each sample:

1. add a pseudovalue of 1e-5 to each non-negative fraction;
2. calculate the centered log-ratio (CLR) of all 53 fractions;
3. estimate feature means and standard deviations in the training fold only;
4. standardize training and held-out samples using the training-fold values;
5. fit a class-balanced ranger probability random forest;
6. use the held-out probability of response as ImmuICBscore.

The fixed model arguments are:

- num.trees = 1200
- mtry = floor(sqrt(53)) = 7
- min.node.size = 8
- class.weights = 0.5 divided by training-fold class prevalence
- probability = TRUE
- five pooled stratified folds

## Reported analysis

The reported pooled analysis contains 334 samples from eight selected ICB
cohorts:

- anti-PD1_SRP070710
- anti-PD1_SRP230414
- anti-PD1_SRP351936
- anti-PD1_ERP105482
- anti-PD1-anti-CTLA4_ERP105482
- anti-PD1_ERP107734
- anti-PD1_ERP117672
- anti-CTLA4-to-anti-PD1_SRP417444

There are 120 responders and 214 non-responders. The archived out-of-fold
predictions give AUC 0.805607, reported as 0.806.

This is a pooled five-fold cross-validation result after exploratory cohort
selection. It is not an independent prospective or treatment-homogeneous
validation estimate. This qualification should accompany the result.

## Reproduce

Required R packages:

~~~r
install.packages(c("data.table", "ranger", "ggplot2"))
~~~

From the repository root:

~~~bash
Rscript analysis/icb/immuicbscore_pooled_fivefold_cv.R \
  analysis/icb/data/icb_immune_fractions.tsv \
  analysis/icb/data/icb_clinical.tsv \
  analysis/icb/reproduced_results \
  8
~~~

The script fixes the fold seed to 20260701 and the forest seeds to 9001 plus
the fold number, matching the selected eight-cohort analysis.

## Input schemas

Fraction table:

- first column: sample identifier;
- remaining columns: the 53 ImmuCellAI 2.0 fractions.

Clinical table:

- Run: sample identifier;
- SRA_study: cohort label;
- ResponseBinary: 1 for responder and 0 for non-responder.

## Archived outputs

The results directory contains the out-of-fold predictions, pooled summary,
and per-cohort AUC table used for the reported result. These files allow direct
verification of the displayed AUC independently of random-forest software
version differences.
