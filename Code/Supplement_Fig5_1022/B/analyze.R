
d <- read_tsv(file.path(data_dir, "B/delta_area.tsv"))
check_num(d$CDF_Area_100bins, analysis$A$areas, "CDF area", 1e-14)
check_num(d$relative_delta_CDF_area, c(NA, diff(analysis$A$areas) / head(analysis$A$areas, -1)), "Relative delta", 1e-14)
analysis$B <- d
