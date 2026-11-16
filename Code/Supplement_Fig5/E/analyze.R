
if (!requireNamespace("survival", quietly = TRUE)) stop("Install survival.")
suppressPackageStartupMessages(library(survival))
d <- context$complete
d$phenotype <- factor(d$phenotype, levels = paste0("IP", 1:4))
d$stage <- factor(d$stage, levels = c("I", "II", "III", "IV"))
d$sex <- factor(d$sex, levels = c("female", "male"))
d$cancer <- factor(d$cancer); d$age10 <- d$age / 10; d$start <- 0
stopifnot(nrow(d) == 4593L, sum(d$OS_event) == 1729L)
archived <- read_tsv(file.path(data_dir, "E/IP4_vs_IP3.tsv"))
arch_models <- read_tsv(file.path(data_dir, "E/model_summary.tsv"))
specs <- list(Unadjusted = character(), Cancer_stratified = "strata(cancer)",
              Clinical_adjusted = c("age10", "sex", "stage", "strata(cancer)"))
contrasts <- comparisons <- diagnostics <- list()
fit_model <- function(terms) coxph(as.formula(paste("Surv(start,OS_time_days,OS_event) ~", paste(terms, collapse = "+"))),
  data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE, singular.ok = FALSE, control = coxph.control(iter.max = 100))
for (id in names(specs)) {
  base <- specs[[id]]
  full <- fit_model(c("phenotype", base)); reduced <- fit_model(if (length(base)) base else "1")
  stopifnot(identical(rownames(full$model), rownames(reduced$model)), full$n == 4593L, full$nevent == 1729L)
  w <- setNames(rep(0, length(coef(full))), names(coef(full)))
  w["phenotypeIP4"] <- 1; w["phenotypeIP3"] <- -1
  beta <- sum(w * coef(full)); se <- sqrt(as.numeric(t(w) %*% vcov(full) %*% w))
  con <- data.frame(model = id, n = full$n, events = full$nevent, HR = exp(beta), CI_low = exp(beta - 1.96 * se),
                    CI_high = exp(beta + 1.96 * se), p = 2 * pnorm(-abs(beta / se)), beta = beta, se = se)
  expected <- archived[archived$model == id, ]
  check_num(unlist(con[4:9]), unlist(expected[4:9]), "Cox contrast", 1e-10)
  df <- length(coef(full)) - length(coef(reduced)); stopifnot(df == 3L)
  lr <- 2 * (tail(full$loglik, 1) - tail(reduced$loglik, 1))
  p <- pchisq(lr, df, lower.tail = FALSE)
  expected <- arch_models[arch_models$model == id, ]
  check_num(lr, expected$IP_LRT_chisq, "Likelihood ratio", 1e-8)
  check_num(p, expected$IP_LRT_p, "Likelihood ratio P", 1e-10)
  contrasts[[id]] <- con; comparisons[[id]] <- data.frame(model = id, chisq = lr, df = df, p = p)
  diagnostics[[id]] <- cox.zph(full)$table
}
analysis$E <- list(data = do.call(rbind, contrasts), comparisons = do.call(rbind, comparisons), PH = diagnostics)
