# Figures 4-6: source implementations and reproduction guide

The three `source_archive/` directories contain **59 full local R scripts**:
47 for Figure 4, five for Figure 5, and seven for Figure 6. They include analysis,
plotting, data preparation, and exploratory variants, not just function stubs.
All 67 R-file instances found under the local Fig4/Fig5/Fig6 folders are covered;
byte-identical duplicates are stored once. Additional distinct workspace
implementations are also included. `source_manifest.json` records original and
archived SHA-256 hashes, normalized source labels, line counts, and imports.

## Which entry points should I use?

- Use the original scripts in `source_archive/` to inspect the historical
  implementation and plotting choices associated with the submitted figures.
- Use the numbered scripts at each module's top level for the configurable
  workflow on tables matching `data/README.md`. These are adaptations, not a
  guarantee of pixel-identical figure reproduction.
- Do not run every archive file in alphabetical order. Several files are
  alternative analyses, interactive notebooks, or exploratory model searches.
  The panel maps in each module README identify the relevant scripts.

The final multi-panel figures were composed from individual plot exports and
illustration assets. Export names can retain old panel/figure numbers. The maps
refer to the supplied final Fig4, Fig5, and Fig6, not those historical names.

## Configure the archived scripts

Personal absolute paths have been replaced with placeholders. Package names and
`hierarchy.mode = "tcell"` use the current ImmuCellAI 2.0 spelling. The numerical
logic of the original scripts has not been rewritten. Output filenames were
renamed consistently; pre-existing historical outputs may need corresponding
filename changes in a separate working copy.

1. Copy `archive_paths.example.R` to `archive_paths.R` and set its absolute roots.
2. From the repository root, run:

   ```bash
   Rscript analysis/prepare_source_archive.R analysis/archive_paths.R
   ```

3. Run the selected script from `analysis/local_sources/<module>/` after placing
   the inputs and prerequisite outputs at the configured paths. For example:

   ```bash
   Rscript analysis/local_sources/age_associated_remodeling/fig5_age_group_immucellai2_abundance.R
   Rscript analysis/local_sources/age_associated_remodeling/fig5_age_group_immucellai2_heatmap.R
   Rscript analysis/local_sources/age_associated_remodeling/fig5_age_group_immucellai2_major_barplot.R
   ```

The preparation command parses every script before writing, refuses to overwrite
an existing destination, and does **not** execute analyses, install packages, or
copy clinical/expression data. Scripts themselves can install a local package
archive when the package is missing, reuse cached fractions, and overwrite their
own result files. Run them in a dedicated copy of the data layout. Historical
interactive scripts may require objects loaded by earlier code blocks.

| Root | Contents |
|---|---|
| `LOCAL_R_ROOT` | `Fig4/`, `Fig5/`, `Fig6/`, `reference_53celltypesTPM20260518.txt`, `MarkerUsedDeconvolution.txt`; optional local R package tarball |
| `LOCAL_CLUSTER_ROOT` | TCGA consensus-class files and matching clinical tables |
| `LOCAL_PROJECT_ROOT` | Original project-level expression and validation inputs used by exploratory scripts |
| `LOCAL_LEGACY_PROJECT_ROOT` | Older server project layout used by archived interactive scripts |
| `LOCAL_LEGACY_AUX_ROOT` | Auxiliary preprocessing inputs used by older interactive scripts |

A cached fraction matrix is not interchangeable with a new run using different
references, marker genes, sample selections, or parameters. Record those choices
and `sessionInfo()` with reproduced results. `manifest.packages` is an inventory
of explicit imports, not a fully resolved dependency lockfile.

## Interpretation and data availability

Figure 4H corresponds to 334 samples, 120 responders and 214 non-responders,
with the locally archived pooled five-fold AUC of 0.8056074766 (rounded to 0.806).
Cohort selection and model exploration occurred before this reported evaluation.
This is not an untouched external validation or a nested evaluation of the full
selection process. The exploratory scripts are retained to make that provenance
visible; re-running a search is not necessary for scoring a fixed cohort set.

Figure 5 has six displayed age groups, with exact interval boundaries documented
in the aging module. Its eight displayed lineage groups are assembled from the
53 states; they are not the model's computational CD4/CD8 hierarchy.

Figure 6's original group comparisons include IQR filtering before some tests;
see the infection README before interpreting P values. Descriptive longitudinal
means and separate status-by-time models answer different questions.

Large expression matrices, patient-level clinical tables, and intermediate model
outputs are not bundled in this code upload. Dataset identifiers and input
schemas are in each module's `data/README.md`; archived scripts additionally
document exact historical filenames. Thus code completeness does not mean that
a fresh clone alone can reproduce every reported number. Figure 4E's background
illustrations are bundled separately under the tumor module's `assets/`.

## Checks

```bash
Rscript analysis/tests/test_source_tools.R
```

This checks syntax, path preparation without analysis execution, age boundaries,
and the complete nonoverlapping 53-state display aggregation. It does not claim
to rerun TCGA deconvolution, cohort searches, or all statistical analyses.
