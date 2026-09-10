# Figure 2: original seven-method benchmark implementations

This directory contains the complete local comparison and plotting scripts,
adapted for portable paths and the current ImmuCellAI 2.0 package name.

## Files and scope

| Script | Figure panel | Cohort |
| --- | --- | --- |
| `compare_gse164522_variable_background_7tools.R` | 2A | Liu, CRC |
| `compare_gse176078_variable_background_7tools.R` | 2B | Wu, BRCA |
| `compare_gse146771_variable_background_7tools.R` | 2C | Zhang, CRC |
| `plot_4datasets_7tools_correlation_pie_matrix_final_ordered.R` | Per-cohort pies | Defaults to the three Figure 2 cohorts |

The original plotting filename is retained. Its fourth-cohort label handling
is also retained, but no fourth dataset is required by the Figure 2 runner.
The plotting script exports individual panel PDFs/PNGs and plot-data tables;
it does not recreate the subsequent manuscript page assembly or manual edits.

Each comparison script includes the full method implementations, input
alignment, cohort-specific state mappings, expected-proportion parsing,
prediction evaluation, cached-result reuse and status logging. These are
archived analysis implementations, not substitute stubs or instructions to
obtain predictions elsewhere. Source-file hashes are recorded in
`source_manifest.tsv`. Imported top-level function bodies are preserved apart
from the package/identifier rename and `tcell` parameter spelling.

One runtime repair recreates the R session's temporary directory after
BayesPrism, whose `run.prism` can remove it. This lets subsequent socket
workers (including CITMIC) start normally without changing numerical updates.
The plotter also accepts archived method identifiers and displays the current
ImmuCellAI 2.0 name.

## What the seven method labels actually run

| Output label | Implementation retained from the local script |
| --- | --- |
| ImmuCellAI 2.0 | `deconvolve_bulk_matrix`, hierarchical VB, `tcell`, UNKNOWN off, 50 iterations, tolerance 1e-6, prior concentrations 10/5/1, pseudo-depth 100000, seed 123 |
| BayesPrism | `new.prism` and `run.prism`; chain length 600, burn-in 500, thinning 2, alpha 1, seed 123, `update.gibbs=FALSE`; first state-level fractions |
| CIBERSORT | Custom linear nu-SVR using `e1071::svm`, nu in 0.25/0.5/0.75, RMSE selection, 1000 ranked reference genes, quantile normalization |
| CITMIC | `CITMIC::CITMIC` with built-in signatures, `log2(bulk+1)`, weighted=TRUE, base=10, damping=0.90; output row normalization |
| DWLS | Local residual-weighted linear regression with nonnegative truncation and normalization; up to six reweighting steps |
| ImmuCellAI | `ImmuCellAI::ImmuCellAI_new` with native signatures and the original GSVA compatibility patch; output row normalization |
| MuSiC | `MuSiC::music.basic` on normalized reference profiles; S=1, Sigma=1e-8, iter.max=1000, nu=1e-4, eps=0.01; weighted-regression fallback on failure |

The CIBERSORT implementation here does **not** invoke CIBERSORTx. A figure label
of CIBERSORTx cannot be established from these scripts. The DWLS implementation
is a custom weighted-regression routine, not an invocation of the official DWLS
package. The MuSiC call uses fixed cell-size and variance inputs rather than
the full multi-donor `music_prop` workflow. These distinctions must accompany
any scientific comparison; this upload preserves historical calculations and
does not validate equivalence to the respective official workflows.

The BayesPrism input is the supplied expression matrix with
`input.type="count.matrix"`, as in the original scripts; there is no additional
TPM-to-count reconstruction. CITMIC and the original ImmuCellAI use their own
signatures. Other methods start with shared, nonzero genes in the supplied
reference; ImmuCellAI 2.0 additionally uses the immune marker panel when at
least 100 panel genes remain. The detailed settings are also tabulated in
`method_settings.tsv`.

## Installation and inputs

Install the current ImmuCellAI 2.0 R package from this repository. Install the
authors' R distributions of BayesPrism, CITMIC, the original ImmuCellAI and
MuSiC in the same R environment. Third-party source packages are not bundled.
The tumor nu-SVR implementation requires `e1071` and `limma`, not the
CIBERSORT/CIBERSORTx package or web service.

```r
install.packages(c("data.table", "e1071", "ggplot2", "RColorBrewer",
                   "magrittr", "stringr", "BiocManager"))
BiocManager::install(c("limma", "GSVA"))
```

`scRNAtoolVis` is optional for the original plotter; it contains a polygon-pie
fallback when that package is unavailable. The original ImmuCellAI adapter
supports both legacy `GSVA::gsva` and the `ssgseaParam` API, and patches its
loaded namespace for the current R process. Installed package files are not
modified. Preserve a separate R session for other analyses.

Input matrices and naming requirements are described in [data/README.md](data/README.md).
Full historical matrices and raw single-cell preprocessing are not included
in this code upload. The packaged atlas is the default; historical reproduction
requires matching reference and marker-file versions as well as mixtures.

## Run the comparison

From the repository root, copy `config.example.R` to `config.R` if custom paths
are needed. Set data/output paths, R library paths, cores and method switches.
All paths in that configuration are relative to the repository root unless
absolute. Defaults select every method and use the packaged atlas/markers.

```sh
# One cohort:
Rscript benchmarks/fig2/run_cohort.R GSE164522 benchmarks/fig2/config.R

# All three cohorts, followed by the original per-cohort plots:
Rscript benchmarks/fig2/run_cohort.R all benchmarks/fig2/config.R

# Replot existing per-cell metric tables:
Rscript benchmarks/fig2/run_cohort.R plot benchmarks/fig2/config.R
```

The config argument is optional. Replace GSE164522 with GSE176078 or GSE146771
for another individual cohort. `sample_limit` and method switches support
small technical checks; a subset run is not a completed seven-method benchmark.
BayesPrism chain settings can be changed for an explicitly labeled smoke test;
the manuscript defaults remain 600/500. Seeds in the original routines are 123.

For direct R use, set `bulk_file`, `reference_file`, `immune_gene_file`,
`out_dir` and `n_cores`, then source the desired comparison script. The runner
additionally records `sessionInfo.txt`, `run_config.rds` and input MD5 checksums,
and reports an error if any enabled method records a failure.

Use a new output directory when changing inputs or settings: the original
scripts reuse existing method fraction files and do not validate their cache
against new inputs. Inspect `run_status.txt` before interpreting results.

## Outputs and interpretation

- Per-method `*_fraction.txt` and `*_metrics.txt` files.
- `all_methods_per_celltype_metrics.txt`, including evaluated state mappings.
- `comparison_summary_7tools.txt`, Pearson/Spearman wide tables and run status.
- In `plots/`: per-cohort PDF/PNG pies, merged plot data and symbol counts.

Mapping functions are embedded in each script as `state_target_mapper`,
`citmic_target_mapper` and `immucellai_target_mapper`. Fractions for mapped
states are summed before evaluation. Pearson/Spearman correlations, RMSE, MAE,
slope and intercept are computed for each target's gradient samples.

An empty mapping is represented by `/`. An undefined correlation with a
nonempty mapping is represented by `x`; in the historical plotter this category
is named `constant_prediction`, but the evaluator can also return NA with
fewer than three finite pairs or constant truth. Interpret it as no correlation
value unless the underlying cause has been checked. The archived summary also
contains `PenalizedMeanPearson`, which replaces undefined correlations with -1;
it is distinct from the mean over valid correlations.

The plotter retains original cell ordering and exclusions (NKT for GSE164522;
NKT, Th1 and Th17 for GSE176078). It does not infer new label mappings or change
method outputs to match the manuscript image.
