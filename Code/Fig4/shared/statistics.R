test_cd <- function(tab, dat, variable, name) {
  res <- tab
  for (i in seq_len(nrow(tab))) {
    dd <- dat[dat$Cell_Type == tab$Cell_Type[i], ]
    if ("Response" %in% names(tab)) dd <- dd[dd$Response == tab$Response[i], ]
    a <- dd$Cell_Ratio[dd[[variable]] == tab$group1[i]]
    b <- dd$Cell_Ratio[dd[[variable]] == tab$group2[i]]
    wt <- wilcox.test(a, b, exact = FALSE)
    res$p_value[i] <- wt$p.value
    res$statistic[i] <- unname(wt$statistic)
    res$n1[i] <- length(a); res$n2[i] <- length(b)
  }
  res$p_adj_BH <- p.adjust(res$p_value, "BH")
  check(paste0(name,"_Wilcoxon_reproduced"), near(res$p_value, tab$p_value) && near(res$statistic, tab$statistic) && near(res$p_adj_BH, tab$p_adj_BH))
  check(paste0(name,"_sample_counts"), identical(res$n1, as.integer(tab$n1)) && identical(res$n2, as.integer(tab$n2)))
}
