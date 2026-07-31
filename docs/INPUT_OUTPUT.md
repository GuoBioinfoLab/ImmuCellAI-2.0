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
`R-package/ImmuneHierDeconv/inst/extdata/reference_53celltypesTPM20260518.txt`.
It is a genes-by-53-states TPM matrix. A custom reference can be supplied, but
its cell-state column names must be compatible with the selected hierarchy.

## Output orientation and interpretation

`state.fraction` is samples by cell states. Fractions are compositional and sum
to one across modeled states. They describe relative immune composition and do
not directly estimate absolute cell counts, tumor purity, or the immune fraction
of a whole tissue.

When UNKNOWN is enabled, UNKNOWN is included in the compositional denominator.
It represents residual expression not explained by the annotated reference and
must not be reported as an immune cell population.

## Recommended quality checks

```r
stopifnot(all(is.finite(fit$state.fraction)))
stopifnot(all(fit$state.fraction >= 0))
stopifnot(max(abs(rowSums(fit$state.fraction) - 1)) < 1e-8)
length(fit$genes.used)
fit$settings
```
