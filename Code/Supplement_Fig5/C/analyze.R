
d <- read_tsv(file.path(data_dir, "C/pac.tsv"))
expected <- read_tsv(file.path(data_dir, "shared/consensus/diagnostics_recomputed.tsv"))
check_table(d, expected, "K", names(d), 1e-14)
stopifnot(d$K[which.min(d[["PAC_0.1_0.9"]])] == 2)
analysis$C <- d
