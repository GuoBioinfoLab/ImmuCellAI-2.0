# Sourced by ../reproduce_Fig4.R in panel order.
# B: the original sample-level analysis is retained; no deduplication,
# cancer adjustment, or clinical filtering is introduced at archival time.
s <- read_tab("B/survival_plot_data.tsv")
check("B_4762_records", nrow(s) == 4762 && !anyDuplicated(s$sample))
check("B_assignments_match_A", identical(s$phenotype, ann$phenotype[match(s$sample, ann$sample)]))
check("B_valid_survival", all(s$status %in% c(0,1)) && all(is.finite(s$time_days) & s$time_days >= 0))
fit <- survfit(Surv(time_days, status) ~ phenotype, data = s)
km <- data.frame(phenotype = sub("phenotype=", "", rep(names(fit$strata), fit$strata)),
                 time_days = fit$time, n_risk = fit$n.risk, n_event = fit$n.event,
                 n_censor = fit$n.censor, survival = fit$surv, lower = fit$lower, upper = fit$upper)
saved_km <- read_tab("B/KM_plot_coordinates.tsv")
check("B_KM_reproduced", identical(km$phenotype, saved_km$phenotype) &&
      all(vapply(names(km)[-1], function(n) near(km[[n]], saved_km[[n]]), logical(1))))
lr <- survdiff(Surv(time_days, status) ~ phenotype, data = s)
lr_result <- data.frame(test = "Log-rank", chisq = lr$chisq, df = length(lr$n)-1,
                        p_value = pchisq(lr$chisq, length(lr$n)-1, lower.tail = FALSE))
check("B_logrank_reproduced", near(lr_result$chisq, read_tab("B/logrank_test.tsv")$chisq))
