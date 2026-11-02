# Sourced by ../reproduce_Fig4.R in panel order.
# C/D: unpaired, two-sided Wilcoxon with exact=FALSE and default continuity
# correction. IQR filtering affects the plot only, never the test population.
cdata <- read_tab("C/violin_plot_data.tsv")
cstats <- read_tab("C/statistics.tsv")
check("C_334_samples_four_cells", nrow(cdata) == 1336 && length(unique(cdata$sample)) == 334)
retained <- rep(FALSE, nrow(cdata))
for (idx in split(seq_len(nrow(cdata)), interaction(cdata$Cell_Type, cdata$Response, cdata$Biopsy_Time))) {
  x <- cdata$Cell_Ratio[idx]
  q <- quantile(x, c(.25, .75), na.rm=TRUE)
  retained[idx] <- x >= q[1]-1.5*diff(q) & x <= q[2]+1.5*diff(q)
}
check("CD_original_plot_filter", identical(retained, cdata$plot_included %in% c(TRUE, "True", "TRUE")))
test_cd(cstats, cdata, "Biopsy_Time", "C")
