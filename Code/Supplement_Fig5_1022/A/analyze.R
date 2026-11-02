
cdf <- read_tsv(file.path(data_dir, "A/cdf_curves.tsv"))
stopifnot(nrow(cdf) == 500L, setequal(cdf$K, 2:6))
areas <- vapply(2:6, function(k) {
  d <- cdf[cdf$K == k, ]
  stopifnot(nrow(d) == 100L, sum(d$pair_count) == 11778231L, all(d$pair_count >= 0))
  expected <- cumsum(d$pair_count) / sum(d$pair_count)
  check_num(d$CDF, expected, "CDF", 1e-14)
  check_num(d$bin_midpoint, (d$bin_left + d$bin_right) / 2)
  sum(expected * (d$bin_right - d$bin_left))
}, numeric(1))
analysis$A <- list(data = cdf, areas = areas)
