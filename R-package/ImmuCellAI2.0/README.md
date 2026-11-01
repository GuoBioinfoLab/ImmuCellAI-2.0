# ImmuCellAI 2.0 R package (v0.2.0)

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
if (!requireNamespace("ImmuCellAI2.0", quietly = TRUE) ||
    getNamespaceVersion("ImmuCellAI2.0") < "0.2.0") {
  stop("Install ImmuCellAI2.0 v0.2.0 or later in a fresh R session; see Installation above.")
}

bulk <- ImmuCellAI2.0::read_expression_matrix("bulk_TPM.txt")

fit <- ImmuCellAI2.0::run_immucellai2(
  bulk = bulk,
  hierarchy.mode = "tcell",
  add.unknown = FALSE,
  n.cores = 8,
  pseudo.depth = 1e5,
  n.iter = 50,
  vb.tol = 1e-4,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  seed = 123
)

ImmuCellAI2.0::write_deconvolution_outputs(fit, "ImmuCellAI2_results")
head(fit$subtype.fraction)
~~~

run_immucellai2() loads the bundled atlas and marker panel, intersects gene
symbols, removes unusable genes, and records all settings in the returned
object. The manuscript default is hierarchy.mode = "tcell" with
add.unknown = FALSE.

## Main outputs

- subtype.fraction: sample-by-53-terminal-subtype matrix used for reporting.
- type.fraction: sample-by-cell-type matrix aggregated at the middle level.
- lineage.fraction: sample-by-lineage matrix aggregated at the top level.
- genes.used: genes retained after intersecting bulk, reference, and markers.
- settings: parameters and dimensions required for reproducibility.
- hierarchy: exact mapping for lineage -> cell type -> subtype.

The fractions are compositional estimates within the modeled immune
compartment. They are not absolute cell counts or whole-tissue fractions.

Exported matrices are `subtype_fraction.txt`, `type_fraction.txt`, and
`lineage_fraction.txt`. In v0.2.0 these replace the old state, subtype, and major
outputs, respectively, without changing their values. In particular, the old
`subtype.fraction` is now `type.fraction`, not the new finest-level result.
See NEWS.md for migration instructions and use a new output directory.
The hierarchy table retains the legacy column names `state_name`, `subtype`,
and `major_lineage` for terminal subtype, parent cell type, and lineage.

## Analysis modes

- tcell: refines CD4/CD8 branches; all other states remain independent.
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

The supported manuscript workflow and `run_immucellai2()` use deterministic
VB inference. The historical Python Gibbs/MCMC implementation is not included
in the current distribution.
