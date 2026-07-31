# ImmuICBscore input schema

The reproduction script expects two processed tables placed in this directory.

## icb_immune_fractions.tsv

- first column: unique sample identifier;
- remaining 53 columns: ImmuCellAI 2.0 relative immune fractions;
- one row per pretreatment ICB sample.

## icb_clinical.tsv

Required columns:

- Run: sample identifier matching the fraction table;
- SRA_study: cohort label;
- ResponseBinary: 1 for responder and 0 for non-responder.

The script selects the eight public SRA/ENA cohorts listed in ../README.md.
Study accession identifiers are provided there so source data can be obtained
under the terms of the original repositories. Processed sample-level clinical
labels are not redistributed here unless their public redistribution status has
been confirmed.
