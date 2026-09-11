# Healthy peripheral-blood aging inputs

The analysis combines the two gene-by-sample TPM matrices `age_bulk.txt` and
`age_bulk2.txt`. Column 1 contains HGNC gene symbols. Sample columns must be
unique across the two files.
The healthy samples used in the compiled metadata originate from public studies
PRJEB14743, PRJEB37238, PRJNA258216, PRJNA509461, PRJNA553703, PRJNA562305,
PRJNA595691, PRJNA662344, PRJNA717024, and PRJNA771014 (SRA studies ERP016409,
ERP120543, SRP045500, SRP173298, SRP214077, SRP219679, SRP241873, SRP281425,
SRP312015, and SRP341241).

Age metadata are read from `age_sample_info.csv` and
`healthy_info2_clean.csv`. Each table must contain a sample identifier and age
in years. Recognized column names are documented in `config.example.R`.

Figure 5 displays six groups: 0-1, 10-20, 20-30, 30-50, 50-70, and 70+.
The exact intervals in the original script are `[0,2)`, `(10,20]`, `(20,30]`,
`(30,50)`, `[50,70)`, and `[70,infinity)`. Ages in `[2,10]`, negative ages,
and missing ages are unassigned. These historical boundaries are retained for
reproducibility, not proposed as a general standard for age binning.
The first metadata table uses `Sample ID` and `Age`; the second uses `Run`
and `Age`. The original scripts retain additional source-batch metadata.
