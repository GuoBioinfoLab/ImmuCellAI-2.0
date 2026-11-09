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
It is a genes-by-53-subtypes TPM matrix. For a custom reference, use
`run_immucellai2_custom(bulk, reference)`. Reference columns can have arbitrary
population names. This workflow uses all shared, informative genes with flat
deterministic Bayesian inference, without marker selection, a biological
hierarchy, or an added unknown component. See the
[custom-reference instructions](API.md#custom-reference-atlas).

## Output orientation and interpretation

For the bundled-reference workflow with UNKNOWN disabled, `subtype.fraction`
has 53 columns by default and 51 with `tissue.type = "non-tumor"`. In the latter
case, TAM and MDSC are absent from the fitted reference and all output levels.
`type.fraction` and `lineage.fraction` are
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

For custom references, `subtype.fraction` instead contains one column for each
supplied reference population, preserving the reference column order. It is not
restricted to 53 populations. No type or lineage summaries are produced, and
the writer exports only `subtype_fraction.txt`. Sample weights describe the
composition relative to the supplied reference, not an absolute cell count.

## Recommended quality checks

```r
stopifnot(all(is.finite(fit$subtype.fraction)))
stopifnot(all(fit$subtype.fraction >= 0))
stopifnot(max(abs(rowSums(fit$subtype.fraction) - 1)) < 1e-8)
length(fit$genes.used)
fit$settings
```
