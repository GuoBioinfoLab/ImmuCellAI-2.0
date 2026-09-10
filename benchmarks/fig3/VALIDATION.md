# Validation record

The portable Figure 3 workflow was checked on 10 September 2026 with R 4.6.0
on Windows.

- All 22 R files under `benchmarks/fig3/` parsed successfully.
- The 18 imported source scripts contained 134 top-level function definitions;
  all 134 were retained after package-name and portable-path adaptation.
- ImmuCellAI 2.0 version 0.1.8 and the six comparator packages were detected.
- The Newman and Monaco preparation stage aligned expression and truth sample
  identifiers successfully.
- Archived numerical result tables were passed through every plotting stage.
  The Xu, Newman, Monaco, runtime, six-cohort summary, and combined A-E outputs
  were generated without errors.
- The six-cohort summary contained 42 expected cohort-by-method combinations.
- A three-sample Newman PBMC smoke test executed all seven methods and produced
  per-cell-type metrics. For this execution-only test, BayesPrism used
  `chain.length = 20` and `burn.in = 10`; manuscript defaults remain
  600 and 500 in the production scripts.
- The smoke test used the package-native LM22 CIBERSORT, native DWLS,
  `MuSiC::music.basic`, native CITMIC, and native ImmuCellAI paths.

The smoke test emitted a MuSiC/snowfall warning about an unrecognized
`--file` command-line option on Windows, but the method completed and its
prediction matrix and metrics were written. Numerical benchmark outputs are
not committed; each full reproduction records `sessionInfo.txt` and
`run_config.rds` in the configured output directory.
