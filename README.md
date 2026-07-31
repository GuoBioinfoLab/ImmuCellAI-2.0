# ImmuCellAI 2.0

ImmuCellAI 2.0 estimates the relative abundance of 53 immune cell states from
bulk RNA-seq expression profiles. The manuscript implementation uses a
reference-based deterministic variational Bayesian (VB) model and targeted
hierarchical refinement of CD4 and CD8 T-cell states.

This repository contains the complete R implementation used for the current
manuscript, the 53-state reference atlas, the 5,510-gene marker panel,
reproducible examples, benchmark utilities, and downstream analysis scripts.
The earlier Python Gibbs-sampling implementation is retained under
[`immucellai2/`](immucellai2/) for compatibility and is explicitly documented
as legacy code.

## Web server

https://guolab.wchscu.cn/ImmuCellAI2/

## Repository layout

| Path | Contents |
|---|---|
| `R-package/ImmuCellAI2.0/` | Complete ImmuCellAI 2.0 R package source |
| `R-package/ImmuCellAI2.0/inst/extdata/` | 53-state reference atlas and 5,510-gene panel |
| `examples/` | Standard run, mode sensitivity, and simulation validation |
| `benchmarks/` | Dataset-independent benchmarking workflow and metadata templates |
| `analysis/` | TCGA immune phenotyping and ImmuICBscore code |
| `docs/` | Algorithm, parameters, input/output, and reproducibility details |
| `immucellai2/` | Legacy Python Gibbs-sampling implementation |

## System requirements

- **Operating system:** Windows, Linux, or macOS. No operating-system-specific
  code is used by the current R implementation.
- **R:** version 4.0.0 or later.
- **Required R packages:** `stats`, `utils`, and `parallel` (included with R;
  all tested at version 4.6.0).
- **Installation helper:** `remotes` is used by the commands below (tested at
  version 2.5.0).
- **Optional packages:** `testthat` 3.3.2 was used for package tests; benchmark
  and downstream analysis scripts list additional packages at the top of each
  script.
- **Hardware:** no non-standard hardware, GPU, or high-performance computing
  environment is required. Multicore CPUs can process samples concurrently.

The package has been tested with R 4.6.0 on Windows 11 x64 and is checked by
GitHub Actions on the current R release using the versioned `windows-latest` and
`ubuntu-latest` hosted environments; exact runner image versions are recorded in
each workflow run.

The recorded demo timing below was obtained on Windows 11 with an Intel Core
Ultra 7 268V CPU and 32 GB RAM; substantially less memory is sufficient for
the bundled four-sample demo.

## Installation

R 4.0 or later is required. From R, install directly from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github(
  "GuoBioinfoLab/ImmuCellAI-2.0",
  subdir = "R-package/ImmuCellAI2.0"
)
```

For a downloaded repository:

```r
remotes::install_local("R-package/ImmuCellAI2.0")
```

Installation from a local source archive took approximately 3 seconds on the
tested desktop. Installation directly from GitHub typically takes 1-3 minutes,
depending mainly on network speed and whether `remotes` is already installed.

Detailed installation and troubleshooting instructions are in
[`docs/INSTALLATION.md`](docs/INSTALLATION.md), and the public API is documented
in [`docs/API.md`](docs/API.md).

## Input

The recommended input is a non-negative TPM matrix with genes in rows and
samples in columns. The first column of a text file must contain gene symbols.

```text
gene    Sample_1    Sample_2
CD3D    42.1        38.7
CD8A    18.2        21.4
MS4A1   12.4        10.8
```

Gene symbols are matched exactly. Duplicate genes are collapsed by their mean.
Counts, FPKM, log-transformed values, and batch-corrected values should not be
mixed in one analysis. See [`docs/INPUT_OUTPUT.md`](docs/INPUT_OUTPUT.md).

## Standard analysis

```r
library(ImmuCellAI2.0)

bulk <- read_expression_matrix("bulk_TPM.txt")

fit <- run_immucellai2(
  bulk = bulk,
  hierarchy.mode = "tcell_only",
  add.unknown = FALSE,
  n.cores = 8,
  seed = 123
)

write_deconvolution_outputs(fit, "ImmuCellAI2_results")
head(fit$state.fraction[, 1:6])
```

`run_immucellai2()` automatically loads the packaged reference and markers,
intersects them with the bulk matrix, and uses the manuscript settings:

```r
immucellai2_defaults()
# hierarchy.mode   = "tcell_only"
# inference.method = "vb"
# add.unknown      = FALSE
# pseudo.depth     = 1e5
# n.iter           = 50
# vb.tol           = 1e-6
# alpha.major      = 10
# alpha.sub        = 5
# alpha.state      = 1
# seed             = 123
```

The full command-line-style example is
[`examples/01_standard_deconvolution.R`](examples/01_standard_deconvolution.R).

## Bundled demo

The repository includes a small synthetic TPM matrix
[`example_data/toy_bulk_TPM.txt`](example_data/toy_bulk_TPM.txt) and its known
mixture proportions
[`example_data/toy_truth_fraction.txt`](example_data/toy_truth_fraction.txt).
After installing the package, run the complete demo from the repository root:

```bash
Rscript examples/01_standard_deconvolution.R \
  example_data/toy_bulk_TPM.txt \
  example_output/toy_ImmuCellAI2 \
  2
```

The demo writes sample-by-cell-state, subtype, and major-lineage fraction
tables plus the genes and settings used. For the bundled input, the primary
state-fraction table contains 4 samples and 53 immune states; all values are
finite and each row sums to 1 when `add.unknown = FALSE`. On the test system
described above, deconvolution took approximately 11 seconds with two CPU
cores. Runtime varies with sample count, overlapping genes, CPU, and storage.

## Output

The returned object contains:

| Field | Orientation | Description |
|---|---|---|
| `state.fraction` | sample x 53 states | Primary reported immune fractions |
| `subtype.fraction` | sample x subtype | Fractions aggregated by internal hierarchy |
| `major.fraction` | sample x major branch | Fractions aggregated by major branch |
| `hierarchy` | table | State-to-subtype-to-major mapping used |
| `genes.used` | vector | Genes retained after all intersections |
| `settings` | list | Exact analysis settings and marker counts |

Rows of `state.fraction` sum to one when `add.unknown = FALSE`. These values are
relative fractions within the modeled immune compartment, not absolute cell
counts or whole-tissue fractions.

## Algorithm summary

For gene `g`, cell state `k`, sample `n`, and iteration `t`, the responsibility
assigned to state `k` is

```math
r_{gnk}^{(t)} =
\frac{\theta_{nk}^{(t-1)} R_{gk}}
{\sum_{j=1}^{K}\theta_{nj}^{(t-1)}R_{gj}}.
```

Expected state counts are

```math
c_{nk}^{(t)} = \sum_g x_{gn} r_{gnk}^{(t)}.
```

The expected counts are combined with hierarchical Dirichlet priors to update
major-lineage, subtype, and state posterior means. Final state fractions are the
product of conditional posterior means along each hierarchy path. Iteration
stops when

```math
\max_k |\theta_{nk}^{(t)}-\theta_{nk}^{(t-1)}| < 10^{-6}
```

or after 50 iterations. Full equations and code links are in
[`docs/ALGORITHM.md`](docs/ALGORITHM.md).

## Modes and UNKNOWN

- `tcell_only`: manuscript default; refines CD4/CD8 branches and leaves other
  immune states as independent branches.
- `flat`: no biologically grouped shrinkage.
- `hierarchical`: broad hierarchy across immune lineages.
- `hybrid`: broader lineage grouping with selected state refinements.
- `add.unknown = TRUE`: adds a residual latent state that competes with the 53
  annotated states. It is not a biological cell type and was not used for the
  principal manuscript results.

Run [`examples/02_compare_four_modes.R`](examples/02_compare_four_modes.R) to
compare `flat` and `tcell_only`, with and without UNKNOWN, on the same data.

## Reproducing validation analyses

The simulation validator expects a table with `Sample`, `TargetCell`, and
`Truth` (or `TargetProportion`) columns. It reports Pearson correlation,
Spearman correlation, RMSE, MAE, bias, slope, and intercept for every target.

```bash
Rscript examples/03_simulated_gradient_validation.R \
  predicted_fractions.tsv truth.tsv validation_results
```

Dataset-specific files are not duplicated when redistribution is restricted or
the file exceeds GitHub's size limit. Required filenames, accession metadata,
and input schemas are documented in [`benchmarks/README.md`](benchmarks/README.md).

## ImmuICBscore reproducibility

The score uses all 53 fractions after centered log-ratio transformation,
training-fold standardization, and a class-balanced `ranger` probability random
forest:

- `num.trees = 1200`
- `mtry = floor(sqrt(53)) = 7`
- `min.node.size = 8`
- class weight = `0.5 / class prevalence`
- five stratified folds; preprocessing is estimated within each training fold

The reported pooled result (`AUC = 0.806`, `N = 334`, 120 responders and 214
non-responders) is a cross-validated estimate across eight selected ICB cohorts,
not an independent prospective validation. See
[`analysis/icb/README.md`](analysis/icb/README.md) and the complete script in
[`analysis/icb/immuicbscore_pooled_fivefold_cv.R`](analysis/icb/immuicbscore_pooled_fivefold_cv.R).

## Reproducibility checklist

Every reported analysis should save:

1. package version and Git commit;
2. reference and marker filenames;
3. number and names of genes used;
4. hierarchy mode and UNKNOWN setting;
5. all VB parameters, seed, and core count;
6. sample mapping and truth-definition files;
7. session information.

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Citation

Please cite the ImmuCellAI 2.0 manuscript when using this implementation. The
final bibliographic record will be added after publication.

## License

MIT License. See [`LICENSE`](LICENSE).
