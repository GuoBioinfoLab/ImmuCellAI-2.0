# Inputs for the original Figure 2 scripts

Provide these gene-by-sample, tab-delimited matrices (first column: gene symbol):

| Panel | Cohort | File |
| --- | --- | --- |
| A | Liu, CRC, GSE164522 | GSE164522_simulate_variable_background.txt |
| B | Wu, BRCA, GSE176078 | GSE176078_simulate_variable_background.txt |
| C | Zhang, CRC, GSE146771 | GSE146771_simulate_variable_background.txt |

These are prepared pseudo-bulk matrices, not GEO raw downloads. Sample names
must be `TargetCell_1` through `TargetCell_40`; the original evaluator derives
the target truth fraction as `step * 0.005`. Renaming arbitrary samples to this
format does not create valid ground truth. Target labels must match the
cohort-specific mappings embedded in the comparison scripts.

The full historical mixture matrices and single-cell preprocessing inputs are
not included here. Recover those files from the original analysis or regenerate
them from the source single-cell data with the documented mixture design.
`benchmarks/create_gradient_pseudobulk.R` generates this input format from
cell-type mean profiles, but exact historical mixtures also require the same
profiles, cell labels, randomization and filtering.

The comparison scripts default to the full 53-state reference and 5,510-gene
panel distributed inside the ImmuCellAI 2.0 R package. Set `reference_file` and
`immune_gene_file` explicitly when reproducing an archived reference version;
record the checksums. Comparator methods use all shared nonzero genes, except
where their own implementation applies an additional signature selection.

Keep local input data and configuration files out of Git. Outputs contain
input checksums and the R session information for reproducibility.
