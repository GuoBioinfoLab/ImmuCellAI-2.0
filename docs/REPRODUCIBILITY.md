# Reproducibility guide

## Record software and source version

```r
packageVersion("ImmuneHierDeconv")
sessionInfo()
system("git rev-parse HEAD", intern = TRUE)
```

## Save the complete settings

```r
fit <- run_immucellai2(bulk, n.cores = 8, seed = 123)
saveRDS(fit$settings, "settings.rds")
writeLines(fit$genes.used, "genes_used.txt")
write_deconvolution_outputs(fit, "results")
capture.output(sessionInfo(), file = "sessionInfo.txt")
```

## Avoid information leakage

For predictive modeling, CLR transformation can be applied before splitting,
but centering and scaling parameters must be estimated using the training fold
only. Feature selection, hyperparameter selection, and cohort selection must
also be nested within the evaluation design or disclosed as exploratory.

## Truth definitions

- Simulated gradients: retain the target-state label and designed proportion
  for every sample.
- Flow cytometry: retain the assay panel, gating definition, unit, denominator,
  and sample mapping.
- Single-cell-derived truth: retain barcode-to-sample mapping, annotation level,
  filtering rules, and whether fractions use all cells or immune cells as the
  denominator.

Do not treat a constant prediction as a missing cell type. In figures, report
constant predictions separately from unavailable/no-reference outputs.

## Suggested repository record

Each analysis directory should include `README.txt`, raw-to-analysis sample
mapping, predictions, per-cell metrics, aggregate metrics, plotting code, and a
session information file.
