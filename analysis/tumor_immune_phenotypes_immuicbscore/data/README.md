# Tumor and immunotherapy input data

Large expression matrices are not stored in GitHub. Configure local paths in
`config.R`.

## TCGA inputs

`TCGA_TPM_symbol.txt` is a gene-by-sample TPM matrix. Column 1 contains HGNC
gene symbols and remaining columns contain TCGA aliquot barcodes. The configured
projects are SKCM, LUAD, LUSC, KIRC, BLCA, HNSC, STAD, ESCA, LIHC, CESC, UCEC,
COAD, READ, and MESO.

`TCGA_clinical.tsv` must contain `case_submitter_id` and `project_id`.
`TCGA_survival.tsv` must contain `sample`, `time_months`, and `status`, where
status is 1 for death and 0 for censoring.

## ICB inputs

The eight treatment-defined cohorts are derived from SRP070710, SRP230414,
SRP351936, ERP105482, ERP107734, ERP117672, and SRP417444. ERP105482 contributes
separate anti-PD-1 and anti-PD-1/anti-CTLA-4 cohorts.

`icb_immune_fractions.tsv` contains sample identifiers in column 1 and 53
ImmuCellAI 2.0 fractions in the remaining columns. `icb_clinical.tsv` contains:

- `Run`: sample identifier;
- `SRA_study`: treatment-defined cohort label;
- `ResponseBinary`: 1 for responder and 0 for non-responder;
- `TreatmentTime`: `Pre` or `On` when longitudinal plots are requested.
