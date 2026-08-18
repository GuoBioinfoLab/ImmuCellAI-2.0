# Benchmark reproduction

## Cohort categories

The manuscript benchmark contains two distinct validation categories.

1. Pseudo-bulk tumor cohorts with known proportions derived from single-cell
   barcodes and source annotations: GSE164522, GSE176078, and GSE146771.
2. Experimentally measured bulk PBMC cohorts with flow-cytometry-derived
   proportions: GSE107019, GSE65133, and GSE107011.

These truth sources must not be described interchangeably. Full provenance is
listed in dataset_metadata.tsv.

## Required files for each cohort

Create one directory per cohort containing:

- bulk_TPM.txt: genes by samples, first column is gene symbol;
- truth_fraction.txt: samples by reported cell types;
- sample_mapping.txt: source sample identifiers and analysis identifiers;
- celltype_mapping.txt: truth labels and prediction columns used for each tool;
- README.txt: filtering, denominator, and provenance notes.

Raw data should be downloaded from the cited public repositories. When a source
file cannot legally be redistributed, record its stable URL, accession,
checksum, and preprocessing command.

## Run ImmuCellAI 2.0

~~~bash
Rscript examples/01_standard_deconvolution.R \
  path/to/bulk_TPM.txt \
  results/ImmuCellAI2 \
  8
~~~

The manuscript settings are tcell VB, UNKNOWN disabled, 50 maximum
iterations, tolerance 1e-6, and Dirichlet concentrations 10, 5, and 1.

## Prediction and truth formats

Prediction table:

~~~text
sample  BGC  Basophil  ...  pDC
S1      0.01 0.001     ... 0.003
S2      0.02 0.002     ... 0.004
~~~

Truth may be wide, as above, or long:

~~~text
Sample  TargetCell  Truth
S1      BGC         0.005
S2      BGC         0.010
~~~

Optional mapping table:

~~~text
TargetCell      PredictionColumns
memory B cell   MBC
Monocytes       cMo+intMo+ncMo
mDCs            monoDC
~~~

The plus sign denotes a sum of prediction columns. An empty
PredictionColumns entry means that a method does not report the target.

## Calculate per-cell metrics

~~~bash
Rscript benchmarks/evaluate_predictions.R \
  results/ImmuCellAI2/state_fraction.txt \
  path/to/truth_fraction.txt \
  results/ImmuCellAI2_evaluation \
  path/to/celltype_mapping.txt
~~~

Outputs:

- aligned_sample_target_values.txt
- per_celltype_metrics.txt
- summary_metrics.txt

The evaluator reports Pearson correlation, Spearman correlation, RMSE, MAE,
bias, slope, and intercept. It also labels each target as valid,
constant_prediction, unavailable, or insufficient_samples.

## Fair comparison rules

- Align sample identifiers before calculating any metric.
- Record the denominator of every truth fraction.
- Do not convert unavailable states to zero.
- Do not calculate a correlation for a constant prediction.
- Distinguish method coverage from prediction accuracy.
- Use the same target mapping for every figure derived from one cohort.
- Report both per-cell metrics and an aggregate summary.
- Record native/default references separately from custom references.
- Record software versions, seeds, core counts, and all non-default arguments.

method_settings.tsv records the manuscript configurations at a concise level.
Dataset-specific state mappings and any deviations must be retained beside the
corresponding analysis.


## Generate a gradient pseudo-bulk benchmark

When cell-type average TPM profiles are available, reproduce the 0.005 to 0.2
target-gradient design with a fixed random background:

~~~bash
Rscript benchmarks/create_gradient_pseudobulk.R \
  celltype_mean_TPM.txt \
  simulation_output \
  target_cells.txt \
  123
~~~

The input is a gene-by-cell-type matrix. The script writes pseudo-bulk TPM,
the complete sample-by-cell truth matrix, the long target truth table, and all
simulation settings. Each target has 40 samples with target proportion
seq(0.005, 0.2, by = 0.005); the remaining mass is randomly distributed over
other states using a fixed seed.

## Reproduce the per-cell correlation figure

Combine per-method metric tables into a tab-delimited file with these columns:

~~~text
Dataset  Method  TargetCell  Pearson  Status
GSE146771  ImmuCellAI 2.0  CD4Tn  0.98  valid
GSE146771  DWLS  CD4Tn  NA  constant_prediction
GSE146771  ImmuCellAI  CD4Tn  NA  unavailable
~~~

Then run:

~~~bash
Rscript benchmarks/plot_per_celltype_correlations.R \
  combined_per_celltype_metrics.tsv \
  Figure_per_celltype_correlations.pdf \
  0.65
~~~

Install the plotting dependency once with:

~~~r
remotes::install_github("junjunlab/scRNAtoolVis")
~~~

The plotting script uses scRNAtoolVis::geom_jjpie. Valid correlations are shown
as Spectral-color pies; constant predictions, unavailable states, and
insufficient sample counts remain visually distinct.
