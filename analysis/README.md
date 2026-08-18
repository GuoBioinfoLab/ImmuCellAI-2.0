# Manuscript case-study analyses

This directory contains the complete downstream analysis and plotting code for
the three ImmuCellAI 2.0 application modules in the manuscript.

| Directory | Case study |
|---|---|
| `tumor_immune_phenotypes_immuicbscore/` | TCGA immune phenotypes, survival, ImmuICBscore, and immunotherapy response |
| `age_associated_remodeling/` | Age-associated immune-cell remodeling in healthy peripheral blood |
| `infectious_disease_states/` | Immune-composition differences in tuberculosis and HIV-associated disease |
| `icb/` | Backward-compatible standalone pooled five-fold ImmuICBscore script and archived outputs |

The large public expression matrices are not duplicated in GitHub. Each module
contains a `data/README.md` with accession identifiers, exact input schemas, and
expected filenames. Copy `config.example.R` to `config.R`, edit the paths, and
run the numbered scripts or `run_all.R` from the repository root.

All deconvolution scripts use the manuscript settings:

```r
hierarchy.mode = "tcell"
add.unknown = FALSE
pseudo.depth = 1e5
n.iter = 50
vb.tol = 1e-6
alpha.major = 10
alpha.sub = 5
alpha.state = 1
seed = 123
```

The case-study scripts are intended to document the complete statistical and
figure-generation workflow. Users may skip the deconvolution step and provide a
sample-by-53-state fraction table conforming to the documented schema.
