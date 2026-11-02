b <- read_tsv(file.path(data_dir, "B/plot_data.tsv"))
check_table(b, analysis$A, "feature")
analysis$B <- b
