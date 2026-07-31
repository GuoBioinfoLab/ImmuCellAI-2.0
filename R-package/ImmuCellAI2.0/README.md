# ImmuCellAI 2.0 R package (v0.1.7)

The `ImmuCellAI2.0` R package implements the ImmuCellAI 2.0 framework
used in the manuscript. It estimates relative fractions of 53 immune cell
states from bulk RNA-seq TPM profiles using deterministic variational Bayesian
inference and targeted hierarchical refinement of CD4 and CD8 T-cell states.

## Installation

From GitHub:

~~~r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github(
  "GuoBioinfoLab/ImmuCellAI-2.0",
  subdir = "R-package/ImmuCellAI2.0"
)
~~~

From a downloaded repository:

~~~r
remotes::install_local("R-package/ImmuCellAI2.0")
~~~

The package imports only base-recommended R packages (stats, utils, and
parallel). The 53-state reference atlas and 5,510-gene marker panel are
bundled in inst/extdata.

## Minimal analysis

The input matrix must contain gene symbols in rows and samples in columns.
Non-negative TPM values are recommended.

~~~r
library(ImmuCellAI2.0)

bulk <- read_expression_matrix("bulk_TPM.txt")

fit <- run_immucellai2(
  bulk = bulk,
  hierarchy.mode = "tcell_only",
  add.unknown = FALSE,
  n.cores = 8,
  pseudo.depth = 1e5,
  n.iter = 50,
  vb.tol = 1e-6,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  seed = 123
)

write_deconvolution_outputs(fit, "ImmuCellAI2_results")
head(fit$state.fraction)
~~~

run_immucellai2() loads the bundled atlas and marker panel, intersects gene
symbols, removes unusable genes, and records all settings in the returned
object. The manuscript default is hierarchy.mode = "tcell_only" with
add.unknown = FALSE.

## Main outputs

- state.fraction: sample-by-53-cell-state matrix used for reporting.
- subtype.fraction: sample-by-subtype aggregated matrix.
- major.fraction: sample-by-major-branch aggregated matrix.
- genes.used: genes retained after intersecting bulk, reference, and markers.
- settings: parameters and dimensions required for reproducibility.
- hierarchy: exact state-to-subtype-to-major mapping.

The fractions are compositional estimates within the modeled immune
compartment. They are not absolute cell counts or whole-tissue fractions.

## Analysis modes

- tcell_only: refines CD4/CD8 branches; all other states remain independent.
- flat: applies no grouped hierarchy and is useful as a sensitivity analysis.
- hierarchical: applies the full broad immune hierarchy.
- hybrid: combines broader grouping with selected state refinements.
- add.unknown = TRUE: adds a latent residual state for unexplained expression;
  the residual is not an annotated biological cell population.

## Complete documentation

See the repository-level documentation:

- docs/INSTALLATION.md: installation, clean-environment test, troubleshooting.
- docs/API.md: complete arguments and returned fields.
- docs/ALGORITHM.md: mathematical formulation and code-to-equation mapping.
- docs/PARAMETERS.md: defaults, interpretation, and sensitivity settings.
- docs/INPUT_OUTPUT.md: matrix formats and validation rules.
- docs/REPRODUCIBILITY.md: reporting checklist.
- examples/: executable standard, sensitivity, and validation workflows.
- benchmarks/: cohort metadata and performance metric scripts.

The older Gibbs/MCMC code paths are retained only for compatibility. The
manuscript workflow and run_immucellai2() use deterministic VB inference.
