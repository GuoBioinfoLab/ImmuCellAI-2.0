# ImmuCellAI 2.0 R package

The `ImmuCellAI2.0` R package performs high-resolution immune-cell deconvolution
of bulk RNA-seq TPM profiles. It estimates relative fractions of 53 immune-cell
subtypes using deterministic variational Bayesian inference and targeted
hierarchical refinement within the CD4 and CD8 T-cell compartments. Output
levels follow lineage -> cell type -> subtype.

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
parallel). The 53-subtype reference atlas and 5,510-gene marker panel are
bundled in inst/extdata.

## Minimal analysis

The input matrix must contain gene symbols in rows and samples in columns.
Non-negative TPM values are recommended.

~~~r
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
add.unknown = FALSE. The default `tissue.type = "tumor"` retains all 53
reference subtypes. For non-tumor samples, use
`run_immucellai2(bulk, tissue.type = "non-tumor")` to fit with 51 reference
subtypes, excluding TAM and MDSC before inference. This retains the standard
marker panel and the requested hierarchy mode. The option only applies to the
bundled reference; the selected setting is recorded in `fit$settings$tissue.type`.

## Main outputs

For the bundled reference:

| R result | Hierarchy level | Exported file |
|---|---|---|
| `fit$subtype.fraction` | Terminal immune subtypes: 53 by default or 51 for non-tumor samples, plus UNKNOWN if enabled | `subtype_fraction.txt` |
| `fit$type.fraction` | Middle level: aggregated cell types | `type_fraction.txt` |
| `fit$lineage.fraction` | Top level: aggregated lineages | `lineage_fraction.txt` |

Additional fields are `fit$genes.used` (retained genes), `fit$settings` (analysis
settings), and `fit$hierarchy` (the mapping used). All three fraction matrices
have samples in rows. Cell-type and lineage values are aggregated subtype
weights, not within-parent conditional proportions.

The fractions are compositional estimates within the modeled immune
compartment. They are not absolute cell counts or whole-tissue fractions.

The hierarchy table uses the column names `state_name`, `subtype`,
and `major_lineage` for terminal subtype, parent cell type, and lineage.

## Custom reference

For users wishing to "retrain" the method for a new dataset or species, such as
livestock immunity studies, prepare a species-matched reference expression
matrix with one column per cell population and a bulk expression matrix with
matching gene identifiers. Use non-negative, non-log-transformed TPM values and
supply both files to `run_immucellai2_custom()` as shown below.

Download the [example bulk matrix](https://raw.githubusercontent.com/GuoBioinfoLab/ImmuCellAI-2.0/main/example_data/custom_reference/bulk_TPM.txt)
and [example reference matrix](https://raw.githubusercontent.com/GuoBioinfoLab/ImmuCellAI-2.0/main/example_data/custom_reference/custom_reference_TPM.txt)
into your R working directory (`getwd()`) before running this example.
The files contain synthetic TPM data for 1,500 genes, three reference populations,
and four bulk samples. The gene and population labels are illustrative, not
biological annotations.

~~~r
fit <- ImmuCellAI2.0::run_immucellai2_custom(
  bulk = "bulk_TPM.txt",
  reference = "custom_reference_TPM.txt"
)
head(fit$subtype.fraction)
ImmuCellAI2.0::write_deconvolution_outputs(fit, "custom_reference_results")
~~~

Both matrices have genes in rows and samples or reference populations in
columns. This workflow uses all shared, informative genes with deterministic
flat Bayesian inference. No marker panel, biological hierarchy, or unknown
component is used. Reference population names are unrestricted. The returned
`subtype.fraction` has one column per reference population; no type or lineage
aggregations are produced. Only `subtype_fraction.txt` is exported.

## Analysis modes

- tcell: refines CD4/CD8 branches; all other subtypes remain independent.
- flat: applies no grouped hierarchy and is useful as a sensitivity analysis.
- add.unknown = TRUE: adds a residual component for unexplained expression;
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
- [Code/](../../Code): manuscript figure analysis and plotting scripts.

The supported manuscript workflow and `run_immucellai2()` use deterministic
VB inference. The historical Python Gibbs/MCMC implementation is not included
in the current distribution.
