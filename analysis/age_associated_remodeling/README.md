# Age-associated immune-cell remodeling

This module contains the healthy peripheral-blood case study and five full
original scripts. See [source preparation](../REPRODUCING_FIGURES.md).

## Figure 5 panel map

| Panel | Original script under `source_archive/` |
|---|---|
| A-C | `fig5_age_group_immucellai2_abundance.R`: merge expression and age metadata, deconvolve 53 states, plot cell-specific boxplots, points and LOESS trends; export unfiltered and plotting-filtered tables. |
| D | `fig5_age_group_immucellai2_major_barplot.R`: aggregate 53 fractions into eight display groups and calculate age-group means. |
| E | `fig5_age_group_immucellai2_heatmap.R`: age-group means, per-state scaling, clipping at -1/1, and ordered heatmap. |

Run the abundance script first, followed by the barplot and heatmap scripts.
They share `<LOCAL_R_ROOT>/Fig5/ImmuCellAI2_age_group_abundance/`.
The other two files, `ABC_plt.r` and `boxplot.r`, are older interactive
preprocessing/plotting implementations with external RData/workbook dependencies;
they are preserved for provenance, not required for the current panels.

The figure contains **six** age groups, not seven. Exact interval boundaries are
listed in [data/README.md](data/README.md). The eight groups in D are T cells,
B cells, ILCs, NK cells, granulocytes, monocytes, macrophages and DCs. The
historical display mapping includes MDSC with monocytes and mast cells with
granulocytes; this is a plotting aggregation, not an ontological claim.

The original A-C script filters plotting outliers within cell-type/age groups
using 1.5 IQR for selected panels, retains CD4Tem/CD8Tem, and applies axis limits
to M0, Bn and CD8Tem. A-B generally retain points except for the CD8Temra
filter. These choices can change plotted summaries and LOESS curves. The
unfiltered fractions, flags, and correlation inputs remain exported. The
numbered adaptation uses different display trimming and summary lines; use the
original script for the historical plotting behavior.

## Configurable workflow

1. `01_prepare_and_deconvolve.R` merges the two TPM matrices, aligns age
   metadata, and runs ImmuCellAI 2.0.
2. `02_plot_age_cell_abundance.R` creates the manuscript abundance panels.
3. `03_plot_age_pattern_heatmap.R` summarizes age-group means, scales each cell
   state across age groups, and draws the ordered heatmap.
4. `04_plot_major_lineage_composition.R` draws the stacked major-lineage plot.

Copy `config.example.R` to `config.R`, edit paths, and run:

```bash
Rscript analysis/age_associated_remodeling/run_all.R \
  analysis/age_associated_remodeling/config.R
```

Display trimming in panel plots does not remove samples from statistical input
or exported tables. The unmodified inferred fractions remain in the
deconvolution output.
