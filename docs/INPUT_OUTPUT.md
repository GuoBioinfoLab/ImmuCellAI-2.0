# Input and output specification

## Bulk expression matrix

- Tab-delimited text or an R numeric matrix.
- Rows are genes; columns are samples.
- For files, the first column contains unique gene symbols.
- Values must be finite and non-negative; TPM is recommended.
- Duplicate genes are collapsed by their arithmetic mean.
- Gene matching is case-sensitive and uses exact symbols.

The standard wrapper intersects bulk genes with both the reference atlas and
the 5,510-gene marker panel. The exact retained genes are returned in
`fit$genes.used` and their count is recorded in `fit$settings`.

## Reference atlas

The packaged file is
`R-package/ImmuCellAI2.0/inst/extdata/reference_53celltypesTPM20260518.txt`.
It is a genes-by-53-subtypes TPM matrix. A custom reference can be supplied, but
its cell-subtype column names must be compatible with the selected hierarchy.

## Output orientation and interpretation

`subtype.fraction` is samples by terminal cell subtypes (53 with the bundled
reference and UNKNOWN disabled). `type.fraction` and `lineage.fraction` are
aggregations at the middle and top levels, respectively. Each matrix contains
sample-level compositional weights, not within-parent conditional proportions.
Fractions sum to one across the modeled populations at each level. They describe
relative immune composition and do
not directly estimate absolute cell counts, tumor purity, or the immune fraction
of a whole tissue.

When UNKNOWN is enabled, UNKNOWN is included in the compositional denominator.
It represents residual expression not explained by the annotated reference and
must not be reported as an immune cell population.

The output names follow lineage -> cell type -> subtype. In the hierarchy table,
`state_name` identifies
the terminal subtype, `subtype` identifies its parent cell type, and
`major_lineage` identifies its lineage. The top-level groups follow the chosen
computational hierarchy; they are not necessarily the nine broad output
categories used in downstream figures.

## Recommended quality checks

```r
stopifnot(all(is.finite(fit$subtype.fraction)))
stopifnot(all(fit$subtype.fraction >= 0))
stopifnot(max(abs(rowSums(fit$subtype.fraction) - 1)) < 1e-8)
length(fit$genes.used)
fit$settings
```
