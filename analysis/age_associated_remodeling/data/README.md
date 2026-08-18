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

The manuscript grouping is: 0-1, 10-20, 20-30, 30-50, 50-70, 70-90, and 90+
years. Samples aged greater than 1 and less than 10 years are unassigned because
that interval was not represented in the reported cohort grouping.
