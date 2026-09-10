# Figure 3 PBMC validation benchmark

This directory contains the complete local R workflow used for the current
manuscript Figure 3: three real PBMC validation cohorts, the seven-method
comparison, the runtime panel, the six-cohort mean-correlation panel, and final
figure assembly.

The imported scripts retain some historical directory and file names such as
`Fig2` and `GSE107011`. Those names were present in the original analysis
workspace and are kept so intermediate files remain traceable. In this README,
the cohorts are identified by their correct accessions.

## Panel-to-code map

| Current panel | Analysis | Principal scripts |
|---|---|---|
| A | Xu PBMC, GSE107019 | `scripts/xu/compare_gse107011_result_tpm_7tools_best_immucellai2.R`, `scripts/xu/compare_gse107019_result_tpm_immucellai2_vs_bayesprism.R`, `scripts/xu/run_gse107011_cibersort_default_lm22.R`, and `scripts/xu/recalculate_plot_gse107011_adjusted_7tools_flatvb.R` |
| B | Alizadeh/Newman PBMC, GSE65133 | `scripts/pbmc/prepare_newman_monaco_inputs.R`, `scripts/pbmc/run_newman_monaco_7tools_compare.R`, and the shared worker |
| C | Larbi/Monaco PBMC, GSE107011 | the same PBMC preparation, runner, and worker scripts |
| D | Runtime efficacy | all scripts under `scripts/runtime/` |
| E | Mean Pearson correlation across six cohorts | `scripts/summary/plot_4datasets_7tools_overall_bar_metrics.R` and `scripts/summary/plot_6datasets_7tools_mean_pearson.R` |
| Full figure | Assemble A-E | `scripts/plots/assemble_current_figure3.R` |

`scripts/plots/make_gse107011_newman_monaco_runtime_combined.R` is the
original A-D assembly script. The current wrapper adds panel E and uses
unambiguous accession-aware display labels.

## What is included

The repository contains the original cohort-specific preprocessing, cell-type
mapping, prediction, evaluation, plotting, and runtime implementations for:

- ImmuCellAI 2.0;
- BayesPrism;
- CIBERSORT;
- CITMIC;
- DWLS;
- ImmuCellAI;
- MuSiC.

This is not a table-only reproduction. The method calls and transformations
used to create each prediction matrix are present in the scripts. See
`method_settings.tsv` before interpreting cross-method results because the
historical workflow used different implementations in different panels. In
particular, PBMC accuracy panels use the package CIBERSORT implementation with
the built-in LM22 signature, whereas the runtime panel uses a local nu-SVR
implementation. The runtime panel therefore measures the archived runtime
configuration, not every possible implementation of each named method.

Private absolute paths and the former package/mode names were replaced with
portable path helpers and the current `ImmuCellAI2.0` / `tcell` names.
These changes do not alter the numerical method definitions. A defensive
wrapper recreates R's temporary directory after BayesPrism because some
BayesPrism versions remove it during cleanup.

## Requirements

Install the current ImmuCellAI 2.0 R package from this repository and the
comparison packages used by the requested stage. A full run requires
`BayesPrism`, `CIBERSORT`, `CITMIC`, `DWLS`, `ImmuCellAI`,
`MuSiC`, `e1071`, `nnls`, `limma`, `GSVA`, `ggplot2`,
`RColorBrewer`, `patchwork`, `magrittr`, and `stringr`.

Use the package release matching the manuscript analysis when exact numerical
reproduction is required. The scripts use the current API name
`hierarchy.mode = "tcell"`.

## Configure

Copy `config.example.R` to `config.R` and edit its paths. The local
`config.R` file is ignored by Git.

The expected input layout is documented in `data/README.md`. Raw controlled
or externally hosted datasets are not duplicated in this repository.

## Run

From the repository root:

```bash
Rscript benchmarks/fig3/run_stage.R prepare_pbmc benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R compare_xu benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R xu_pairwise benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R evaluate_xu benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R compare_pbmc benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R runtime benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R runtime_summary benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R plot_xu benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R plot_pbmc benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R summarize benchmarks/fig3/config.R
Rscript benchmarks/fig3/run_stage.R compose benchmarks/fig3/config.R
```

The entire workflow can be launched with:

```bash
Rscript benchmarks/fig3/run_stage.R all benchmarks/fig3/config.R
```

A full BayesPrism and seven-method run can take substantial time. To regenerate
plots from completed numerical results, use the `plots` stage. The
`plot_runtime` stage accepts an existing combined timing table; the
`runtime_summary` stage instead requires both freshly generated scenario
summaries.

## Reproducibility notes

- The Xu comparison evaluates four ImmuCellAI 2.0 configurations and the final
  displayed result uses flat VB with `add.unknown = FALSE`.
- The Newman workflow uses flat VB with `add.unknown = FALSE`.
- The Monaco workflow preserves the historical setting: T-cell hierarchical VB
  with `add.unknown = TRUE`.
- BayesPrism uses `chain.length = 600`, `burn.in = 500`, thinning 2,
  alpha 1, and seed 123 in the full benchmark.
- Runtime mixtures use seed 20260604 and Dirichlet shape 0.7 for each of the 53
  states; ImmuCellAI 2.0 uses T-cell hierarchical VB, no UNKNOWN state,
  concentrations 10/5/1, tolerance 1e-6, and seed 123.
- The runtime plot never silently invents timings. Archived constants are
  available only when `use_archived_timing = TRUE`; the default is to stop
  if no measured summary table exists.
- Figure 2 tumor-cohort metrics must be generated by
  `benchmarks/fig2/` and supplied through `fig2_results_dir`.

See `VALIDATION.md` for the checks performed on the uploaded workflow and
`source_manifest.tsv` for SHA-256 hashes of the local source scripts before
portable-path adaptations.
