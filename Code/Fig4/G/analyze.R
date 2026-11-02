# Sourced by ../reproduce_Fig4.R in panel order.
# G: phenotype composition; F provides the frozen per-sample quartile labels.
gdata <- read_tab("G/quartile_plot_data.tsv")
counts <- table(quartiles, fdata$cluster)
for (i in seq_len(nrow(gdata))) {
  k <- gdata$score_quartile[i]; cl <- gdata$cluster[i]
  check(paste0("G_",i), counts[k,cl] == gdata$N[i] &&
        near(counts[k,cl]/sum(counts[k,]), gdata$Proportion[i]))
}
