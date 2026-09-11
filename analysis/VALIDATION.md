# Source-upload checks (2026-09-11)

- All 67 R-file instances under the original local Fig4/Fig5/Fig6 directories
  are represented by source hashes in `source_manifest.json`; byte-identical
  duplicates are not stored twice. The manifest contains 59 implementations,
  including two additional workspace scripts. All 59 original hashes and all
  59 normalized archive hashes were checked.
- All 85 R scripts under `analysis/` parsed successfully.
- `Rscript analysis/tests/test_source_tools.R` passed: six-group interval
  boundaries, nonoverlapping aggregation of all 53 states into eight display
  groups, path substitution for all 59 archived scripts, refusal to overwrite
  prepared sources, and preparation without executing an analysis.
- The local historical Figure 4H prediction table was checked independently
  using the rank-based AUC formula: N=334, R=120, NR=214,
  AUC=0.805607476635514. This checks stored predictions, not model retraining.
- The local historical Figure 5 fraction matrix contained 327 samples and 53
  states. The shared six-group function matched all historical metadata labels.
  Counts were 7, 43, 124, 85, 60, and 8 for 0-1, 10-20, 20-30, 30-50,
  50-70, and 70+, respectively.
- Source files were checked for personal absolute path prefixes and common
  credential patterns before upload. Raw clinical and expression tables were
  not included in the upload.

These are source/provenance and focused regression checks, not an end-to-end
reproduction of all figures. Full deconvolution, model searches, Figure 6
statistical analyses, and complete visual regeneration were not rerun for this
upload. The full plotting/model dependency stack was not available in the
verification R session. Follow each module's input/dependency documentation to
run those analyses in the intended environment.
