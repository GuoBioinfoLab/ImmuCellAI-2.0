# Infectious-disease case-study inputs

`exp.txt` is the shared gene-by-sample TPM matrix. Column 1 contains HGNC gene
symbols and the other columns use public run accessions as sample identifiers.

`PRJNA_combined.csv` supplies tuberculosis metadata and must contain `Run` and
`tb_status`, with `Active TB` and `Latent TB` as the comparison labels.The tuberculosis samples are from PRJNA638653 and PRJNA395234.

`PRJNA683803_combined.csv` supplies HIV-associated clinical labels for PRJNA683803 and must
contain `Run`, `outcome`, and `treatment`. `PRJNA683803.xlsx` supplies the
longitudinal fields `Run`, `Patient_ID`, `grouping`, and `outcome`. Grouping
suffixes `_0`, `_1`, `_4`/`_4i`, and `_8` are mapped to D0, D1, D4, and D8.

The source publications and public repository records should be used to obtain
the expression and clinical tables. No patient-level restricted data are
distributed in this repository.
