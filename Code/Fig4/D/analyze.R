# Sourced by ../reproduce_Fig4.R in panel order.
# D: on-treatment response comparison; C provides the parent sample table.
ddata <- read_tab("D/violin_plot_data.tsv")
dstats <- read_tab("D/statistics.tsv")
check("D_90_samples_four_cells", nrow(ddata) == 360 && length(unique(ddata$sample)) == 90)
check("CD_subset", identical(ddata, cdata[cdata$Biopsy_Time == "On", , drop=FALSE]) ||
      isTRUE(all.equal(ddata, cdata[cdata$Biopsy_Time == "On", , drop=FALSE], check.attributes=FALSE)))
test_cd(dstats, ddata, "Response", "D")
