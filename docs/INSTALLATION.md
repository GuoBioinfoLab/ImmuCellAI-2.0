# Installation and environment

## Supported environment

- R 4.0 or later
- Windows, Linux, or macOS
- Base package dependencies: stats, utils, and parallel
- Optional test dependency: testthat

The manuscript implementation is the R package under
R-package/ImmuCellAI2.0. The historical Python code under immucellai2/
does not implement the current deterministic variational Bayesian algorithm.

## Install from GitHub

~~~r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "GuoBioinfoLab/ImmuCellAI-2.0",
  subdir = "R-package/ImmuCellAI2.0",
  upgrade = "never"
)
~~~

Restart R and verify the installation:

~~~r
library(ImmuCellAI2.0)
packageVersion("ImmuCellAI2.0")
immucellai2_defaults()
dim(load_immucellai2_reference())
length(load_immucellai2_markers())
~~~

The final two commands should report 53 reference columns and 5,510 marker
genes.

## Install from a downloaded repository

Set the working directory to the repository root and run:

~~~r
remotes::install_local("R-package/ImmuCellAI2.0", upgrade = "never")
~~~

On Windows, installing this package does not require compilation because it
contains only R code. Rtools is therefore not required for ordinary use.

## Clean-environment verification

~~~r
dir.create("temporary_R_library", showWarnings = FALSE)
.libPaths(c(normalizePath("temporary_R_library"), .libPaths()))

remotes::install_local("R-package/ImmuCellAI2.0", lib = .libPaths()[1], dependencies = FALSE, upgrade = "never")

library(ImmuCellAI2.0, lib.loc = .libPaths()[1])
stopifnot(ncol(load_immucellai2_reference()) == 53L)
stopifnot(length(load_immucellai2_markers()) == 5510L)
~~~

## Run the included smoke example

From the repository root:

~~~bash
Rscript examples/00_make_toy_bulk.R
Rscript examples/01_standard_deconvolution.R \
  example_data/toy_bulk_TPM.txt \
  example_results \
  2
~~~

The output directory should contain:

- state_fraction.txt
- subtype_fraction.txt
- major_lineage_fraction.txt
- hierarchy_used.txt
- settings.rds
- genes_used.txt
- sessionInfo.txt

## Common installation problems

### Package not found after installation

The package source has not been installed into the active library. Re-run the
installation command and inspect .libPaths().

### GitHub download fails behind a proxy

Download the repository as a ZIP archive from GitHub and use the local
installation command. The package does not download reference files at runtime;
the reference atlas and marker panel are bundled.

### The package installs but no genes overlap

Use current HGNC gene symbols in the first column of the input file. Ensembl
identifiers must be mapped to symbols before analysis. Gene matching is exact
and case-sensitive.

### A large analysis appears to use only one core

n.cores parallelizes independent samples. The VB updates within one sample are
sequential, so a single-sample run does not benefit substantially from a large
core count.
