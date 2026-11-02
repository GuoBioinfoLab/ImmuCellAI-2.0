
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
counts <- read_tsv(file.path(data_dir, "D/cancer_IP_counts.tsv"))
fresh <- aggregate(cbind(n = rep(1L, nrow(context$patients)), events = context$patients$OS_event),
                   context$patients[c("cancer", "phenotype")], sum)
check_table(counts, fresh, c("cancer", "phenotype"))
d <- read_tsv(file.path(data_dir, "D/plot_percentages.tsv"))
check_num(d$percent, 100 * counts$n / ave(counts$n, counts$cancer, FUN = sum), "Cancer percentages")
test <- chisq.test(xtabs(n ~ cancer + phenotype, counts), correct = FALSE)
reference <- jsonlite::fromJSON(file.path(data_dir, "D/distribution_test.json"))
value <- sqrt(unname(test$statistic) / (4746 * 3))
check_num(c(unname(test$statistic), value), c(reference$chi_square, reference$cramers_V), "Distribution test", 5e-5)
analysis$D <- list(data = d, test = test, cramers_V = value)
