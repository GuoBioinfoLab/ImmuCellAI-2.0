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
It is a genes-by-53-states TPM matrix. A custom reference can be supplied, but
its cell-state column names must be compatible with the selected hierarchy.

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

The output names follow lineage -> cell type -> subtype. Hierarchy-table column
names are retained for compatibility with custom inputs: `state_name` identifies
the terminal subtype, `subtype` identifies its parent cell type, and
`major_lineage` identifies its lineage. The top-level groups follow the chosen
computational hierarchy; they are not necessarily the nine broad output
categories used in downstream figures. No group assignments changed in v0.2.0.

## Upgrading from v0.1.9 or earlier

| Old result field / file | New result field / file | Level |
|---|---|---|
| `state.fraction` / `state_fraction.txt` | `subtype.fraction` / `subtype_fraction.txt` | Terminal subtype |
| `subtype.fraction` / `subtype_fraction.txt` | `type.fraction` / `type_fraction.txt` | Cell type |
| `major.fraction` / `major_lineage_fraction.txt` | `lineage.fraction` / `lineage_fraction.txt` | Lineage |

This is a naming change only. Old aliases are not returned because the old and
new meanings of `subtype.fraction` conflict. To migrate a previously saved fit
without rerunning inference, rename its fields in one operation:

```r
stopifnot(all(c("state.fraction", "subtype.fraction", "major.fraction") %in% names(fit)))
old_names <- c("state.fraction", "subtype.fraction", "major.fraction")
names(fit)[match(old_names, names(fit))] <-
  c("subtype.fraction", "type.fraction", "lineage.fraction")
```

Use a new output directory. The writer rejects legacy result fields and
directories containing `state_fraction.txt` or `major_lineage_fraction.txt` to
avoid overwriting a middle-level `subtype_fraction.txt` with terminal results.
Historical `source_archive/` scripts and recorded data retain their original
names; their recorded package version should be used for exact reproduction.

## Recommended quality checks

```r
stopifnot(all(is.finite(fit$subtype.fraction)))
stopifnot(all(fit$subtype.fraction >= 0))
stopifnot(max(abs(rowSums(fit$subtype.fraction) - 1)) < 1e-8)
length(fit$genes.used)
fit$settings
```
