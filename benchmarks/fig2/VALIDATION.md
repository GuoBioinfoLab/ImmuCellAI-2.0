# Import validation (2026-09-10)

## Checks completed

- All six R files in this directory parse successfully.
- Each of the three comparison scripts preserves all 24 original top-level
  function definitions. Parsed function expressions were compared after only
  package/identifier and hierarchy-mode spelling substitutions.
- A technical run used the first three samples of the existing GSE146771
  variable-background simulation, one core, and a fresh output directory.
  All seven enabled methods completed and wrote prediction/metric files.
- For that technical run only, BayesPrism used 20 iterations with a burn-in
  of 10. Uploaded defaults remain 600/500. This check is not a recomputation
  of the published correlations or a scientific accuracy benchmark.
- The original plotting script processed archived metric tables for all three
  cohorts. All seven method labels were retained for each cohort, and six
  nonempty per-cohort PDF files (high/low variants) plus PNGs were generated.
  The polygon-pie fallback was exercised because scRNAtoolVis was unavailable.
- Git whitespace/error checks passed. Raw inputs, generated outputs, the local
  test R library and machine-specific configurations are excluded from Git.

## Tested environment

Windows, R 4.6.0 (ucrt). ImmuCellAI 2.0 version 0.1.8 was installed from this
repository into an isolated library for the final seven-method run.

| Package | Version |
| --- | --- |
| ImmuCellAI2.0 | 0.1.8 |
| BayesPrism | 2.2.3 |
| CITMIC | 0.1.3 |
| ImmuCellAI | 0.1.0 |
| MuSiC | 1.0.0 |
| e1071 | 1.7-17 |
| limma | 3.68.4 |
| GSVA | 2.6.2 |
| ggplot2 | 4.0.3 |

These are the import-test versions, not a claim about the historical versions
used to generate every manuscript panel. The current BayesPrism run removed
the session temporary directory; restoring it before later methods resolved
the observed CITMIC socket-worker startup failure.

The scripts retain historical custom and reduced comparator implementations,
including MuSiC's regression fallback. Successful execution does not establish
their equivalence to full official workflows. See README.md and
method_settings.tsv for the exact implementations and interpretation limits.
