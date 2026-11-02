# ImmuCellAI 2.0

ImmuCellAI 2.0 is a reference-based method for high-resolution immune-cell
deconvolution of bulk RNA-seq TPM profiles. It estimates the relative abundance
of 53 immune-cell subtypes using a curated reference atlas and deterministic
variational Bayesian (VB) inference, with targeted hierarchical refinement
within the CD4 and CD8 T-cell compartments.

Results follow the hierarchy **lineage -> cell type -> subtype**. The primary
53-subtype result is `fit$subtype.fraction`; aggregated cell-type and lineage
results are available as `fit$type.fraction` and `fit$lineage.fraction`.

This repository contains the complete R implementation used for the current
manuscript, the 53-subtype reference atlas, the 5,510-gene marker panel,
reproducible examples, benchmark utilities, and downstream analysis scripts.

## Web server

https://guolab.wchscu.cn/ImmuCellAI2/

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

R 4.0 or later is required. The package imports only the R packages `stats`,
`utils`, and `parallel`, which are included with R.

### Single-file download for offline or intranet servers

Users who do not want the complete repository can download only the current R
source package:

**[Download ImmuCellAI2.0_0.3.0.tar.gz](https://github.com/GuoBioinfoLab/ImmuCellAI-2.0/releases/download/v0.3.0/ImmuCellAI2.0_0.3.0.tar.gz)**

All published package versions are available on the
**[GitHub Releases page](https://github.com/GuoBioinfoLab/ImmuCellAI-2.0/releases/latest)**.
The file named `ImmuCellAI2.0_<version>.tar.gz` is the complete installable R
package, including the R code, 53-subtype reference atlas, and marker-gene panel.
Cloning or downloading the rest of the repository is not required.

After transferring this single file to an intranet server, install it from R:

```r
install.packages(
  "/path/to/ImmuCellAI2.0_0.3.0.tar.gz",
  repos = NULL,
  type = "source"
)
library(ImmuCellAI2.0)
```

Alternatively, install it from a shell:

```bash
R CMD INSTALL ImmuCellAI2.0_0.3.0.tar.gz
```

### Install directly from the repository

On an internet-connected computer, the package can also be installed directly
from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github(
  "GuoBioinfoLab/ImmuCellAI-2.0",
  subdir = "R-package/ImmuCellAI2.0"
)
```

For a downloaded complete repository:

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
bulk <- ImmuCellAI2.0::read_expression_matrix("bulk_TPM.txt")

fit <- ImmuCellAI2.0::run_immucellai2(
  bulk = bulk,
  hierarchy.mode = "tcell",
  add.unknown = FALSE,
  n.cores = 8,
  seed = 123
)

ImmuCellAI2.0::write_deconvolution_outputs(fit, "ImmuCellAI2_results")
head(fit$subtype.fraction[, 1:6])
```

`run_immucellai2()` automatically loads the packaged reference and markers,
intersects them with the bulk matrix, and uses the following default settings:

```r
ImmuCellAI2.0::immucellai2_defaults()
# hierarchy.mode   = "tcell"
# inference.method = "vb"
# add.unknown      = FALSE
# pseudo.depth     = 1e5
# n.iter           = 50
# vb.tol           = 1e-4
# alpha.major      = 10
# alpha.sub        = 5
# alpha.state      = 1
# seed             = 123
```

The three settings below control population grouping, inference, and the
optional unknown component in the bundled-reference workflow:

| Setting | Values | Meaning |
|---|---|---|
| `hierarchy.mode` | `"tcell"` (default), `"flat"` | `"tcell"` applies an internal hierarchy within the CD4 and CD8 T-cell branches. `"flat"` places known populations in one common pool without biological hierarchy grouping. |
| `inference.method` | `"vb"` (standard workflow); `"mcmc"` in the lower-level API | `"vb"` uses deterministic posterior-mean updates with the selected hierarchy mode, including `"flat"`. `"mcmc"` is the legacy stochastic Gibbs-sampling method. |
| `add.unknown` | `FALSE` (default), `TRUE` | `FALSE` estimates only the supplied known populations. `TRUE` jointly estimates an additional `UNKNOWN` residual weight; known and unknown weights together sum to one. Set `TRUE` only when the sample is known to contain components not represented in the reference. |

`run_immucellai2()` fixes the inference method to `"vb"`; `inference.method`
is reported in its settings but is **not an argument to this function**.
Alternative inference methods are available through `deconvolve_bulk_matrix()`;
see [`docs/API.md`](docs/API.md#inference-arguments).

The full command-line-style example is
[`examples/01_standard_deconvolution.R`](examples/01_standard_deconvolution.R).

## Example workflow

The repository includes a small synthetic TPM matrix
[`example_data/toy_bulk_TPM.txt`](example_data/toy_bulk_TPM.txt) and its known
mixture proportions
[`example_data/toy_truth_fraction.txt`](example_data/toy_truth_fraction.txt).
After installing the package and downloading the repository, run the following
code directly in the R or RStudio console. Replace the path in `setwd()` with
your repository folder, which contains `example_data/` and `examples/`.
Run the complete block below. The `ImmuCellAI2.0::` prefix calls functions from
the installed package without requiring a separate `library()` call.

```r
setwd("path/to/ImmuCellAI-2.0")

bulk <- ImmuCellAI2.0::read_expression_matrix("example_data/toy_bulk_TPM.txt")

fit <- ImmuCellAI2.0::run_immucellai2(
  bulk = bulk,
  hierarchy.mode = "tcell",
  add.unknown = FALSE,
  n.cores = 2
)

dim(fit$subtype.fraction)
# [1] 4 53

dim(fit$type.fraction)
dim(fit$lineage.fraction)

fit$subtype.fraction[, 1:6, drop = FALSE]
rowSums(fit$subtype.fraction)

output_dir <- "example_output/toy_ImmuCellAI2"
ImmuCellAI2.0::write_deconvolution_outputs(fit, output_dir)
saveRDS(fit$settings, file.path(output_dir, "settings.rds"))
writeLines(fit$genes.used, file.path(output_dir, "genes_used.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
```

The console displays the dimensions at all three hierarchy levels, the first
six subtype fractions for each sample, and the subtype row sums (approximately
1). The complete fraction tables
and run settings are saved in `example_output/toy_ImmuCellAI2/`. The known mixture
proportions are provided for comparison and are not used as inference input.
On the test system described above, the four-sample demo took approximately
11 seconds with two CPU cores; runtime varies across systems.

## Output

For the bundled reference, the three fraction matrices have samples in rows and the corresponding
populations in columns:

| R result | Hierarchy level | Exported file |
|---|---|---|
| `fit$subtype.fraction` | Finest level: 53 terminal immune subtypes | `subtype_fraction.txt` |
| `fit$type.fraction` | Middle level: aggregated cell types | `type_fraction.txt` |
| `fit$lineage.fraction` | Top level: aggregated lineages | `lineage_fraction.txt` |

The object also includes `fit$hierarchy` (the mapping used), `fit$genes.used`
(retained genes), and `fit$settings` (analysis settings and marker counts).

Rows of `subtype.fraction` sum to one when `add.unknown = FALSE`. These values are
relative fractions within the modeled immune compartment, not absolute cell
counts or whole-tissue fractions.

`write_deconvolution_outputs()` writes the three fraction files listed above
and `hierarchy_used.txt`. Cell-type and lineage results are sums of the terminal
subtype weights, not within-parent conditional proportions. Their columns follow
the selected computational hierarchy; lineage groups are not necessarily the
nine broad categories used for downstream reporting.

## Custom reference

For users wishing to "retrain" the method for a new dataset or species, such as
livestock immunity studies, prepare a species-matched reference expression
matrix with one column per cell population and a bulk expression matrix with
matching gene identifiers. Use non-negative, non-log-transformed TPM values and
supply both files to `run_immucellai2_custom()` as shown below.

Supply a bulk matrix and a reference matrix with genes in rows and samples or
reference populations in columns. The first column of each tab-delimited file
contains matching gene identifiers. Non-negative TPM values on a compatible scale are
recommended; reference population names can be chosen freely.

Example files: [bulk TPM matrix](example_data/custom_reference/bulk_TPM.txt) and
[custom reference TPM matrix](example_data/custom_reference/custom_reference_TPM.txt).
This synthetic demo contains 1,500 genes, three reference populations, and four
bulk samples. Synthetic gene and population names are for
demonstration only. Run the following from the repository root:

```r
fit <- ImmuCellAI2.0::run_immucellai2_custom(
  bulk = "example_data/custom_reference/bulk_TPM.txt",
  reference = "example_data/custom_reference/custom_reference_TPM.txt"
)
head(fit$subtype.fraction)
ImmuCellAI2.0::write_deconvolution_outputs(fit, "custom_reference_results")
```

This workflow uses all shared genes with positive expression in both inputs. It applies deterministic flat Bayesian posterior-mean inference.
The result contains one column per supplied reference population, with
sample-specific mixture weights summing to one. Only `subtype_fraction.txt` is
exported; no cell-type or lineage aggregations are inferred. The genes used and
run settings are available in `fit$genes.used` and `fit$settings`.

See [custom-reference details](docs/API.md#custom-reference-atlas) and the
[example script](examples/04_custom_reference.R).

## Citation

Please cite the ImmuCellAI 2.0 manuscript when using this implementation. The
final bibliographic record will be added after publication.

## License

MIT License. See [`LICENSE`](LICENSE).
