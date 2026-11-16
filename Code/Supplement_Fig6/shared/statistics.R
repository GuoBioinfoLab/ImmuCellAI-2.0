read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                                     fileEncoding = "UTF-8-BOM", quote = "\"", comment.char = "")
table_key <- function(x, cols) do.call(paste, c(x[cols], sep = "|"))
check_num <- function(a, b, label = "Numeric comparison", tolerance = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) stop(label, ": missingness or length differs")
  finite <- is.finite(a) & is.finite(b)
  if (any(is.infinite(a) != is.infinite(b)) || any(a[!finite & !is.na(a)] != b[!finite & !is.na(b)]))
    stop(label, ": nonfinite values differ")
  if (any(abs(a[finite] - b[finite]) > tolerance)) stop(label, ": max difference ", max(abs(a[finite] - b[finite])))
  invisible(TRUE)
}
check_table <- function(a, b, keys, cols = intersect(names(a), names(b)), tolerance = 1e-12) {
  ka <- table_key(a, keys); kb <- table_key(b, keys)
  stopifnot(!anyDuplicated(ka), !anyDuplicated(kb), setequal(ka, kb))
  b <- b[match(ka, kb), , drop = FALSE]
  for (column in cols) {
    if (is.numeric(a[[column]]) && is.numeric(b[[column]])) check_num(a[[column]], b[[column]], column, tolerance)
    else stopifnot(identical(as.character(a[[column]]), as.character(b[[column]])))
  }
  invisible(TRUE)
}
auc_rank <- function(y, score) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  stopifnot(n1 > 0, n0 > 0, all(is.finite(score)), all(y %in% c(0, 1)))
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

check_equal <- function(actual, expected, columns = intersect(names(actual), names(expected)), tolerance = 1e-8) {
  stopifnot(nrow(actual) == nrow(expected))
  for (column in columns) {
    a <- actual[[column]]; b <- expected[[column]]
    if (is.numeric(a) && is.numeric(b)) {
      stopifnot(identical(is.na(a), is.na(b)))
      keep <- !is.na(a); stopifnot(all(abs(a[keep] - b[keep]) <= tolerance * pmax(1, abs(b[keep]))))
    } else stopifnot(identical(as.character(a), as.character(b)))
  }
}
prepare_patients <- function(d) {
  stopifnot(nrow(d) == 4593L, !anyDuplicated(d$patient), sum(d$OS_event) == 1729L,
            length(unique(d$cancer)) == 14L, all(d$OS_time_days > 0))
  d$sex <- factor(d$sex, levels = c("female", "male"))
  d$stage <- factor(d$stage, levels = c("I", "II", "III", "IV"))
  d$cancer <- factor(d$cancer); d$age10 <- d$age / 10; d$months <- d$OS_time_days / 30.4375
  stopifnot(!anyNA(d[, c("age10", "sex", "stage", "cancer")]))
  d
}
prepare_analysis <- function(data_dir) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Install survival.")
  suppressPackageStartupMessages(library(survival))
  list(patients = prepare_patients(read_tsv(file.path(data_dir, "shared/survival/patient_scores_and_groups.tsv"))),
       verification = list(patients = 4593L, events = 1729L, RF_refit = FALSE))
}
analyze_survival <- function() {
  tab <- function(p) read_tsv(file.path(data_dir, p))
d <- prepare_patients(tab('shared/survival/patient_scores_and_groups.tsv'))
specs <- data.frame(panel=strsplit('CDEFGH','')[[1]],
  split=rep(c('pooled_median','within_cancer_median'),each=3),
  marker=rep(c('ImmuICBscore','Major_RF','CD45'),2),
  column=rep(c('ImmuICBscore','Major_RF_score','CD45_TPM'),2))
results <- list(); ph_rows <- list(); curves <- risks <- list(); warnings <- character()
for(i in seq_len(nrow(specs))) {
  s <- specs[i,]
  cutoff <- if(s$split=='pooled_median') rep(median(d[[s$column]]),nrow(d)) else ave(d[[s$column]],d$cancer,FUN=median)
  d$group <- factor(ifelse(d[[s$column]]>cutoff,'High','Low'),levels=c('Low','High'))
  stopifnot(identical(as.character(d$group),d[[paste0(s$marker,'_',s$split,'_group')]]))
  lr <- survdiff(Surv(OS_time_days,OS_event)~group,data=d,rho=0)
  slr <- survdiff(Surv(OS_time_days,OS_event)~group+strata(cancer),data=d,rho=0)
  raw <- coxph(Surv(OS_time_days,OS_event)~group,data=d,ties='efron',x=TRUE)
  adjusted <- withCallingHandlers(coxph(Surv(OS_time_days,OS_event)~group+age10+sex+stage+strata(cancer),
      data=d,ties='efron',x=TRUE,singular.ok=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  ss <- summary(adjusted); s0 <- summary(raw); ph <- cox.zph(adjusted)$table
  result <- data.frame(split=s$split,marker=s$marker,n=nrow(d),events=sum(d$OS_event),
    low_n=sum(d$group=='Low'),high_n=sum(d$group=='High'),
    low_events=sum(d$OS_event[d$group=='Low']),high_events=sum(d$OS_event[d$group=='High']),
    logrank_chisq=lr$chisq,logrank_p=pchisq(lr$chisq,1,lower.tail=FALSE),
    cancer_stratified_logrank_chisq=slr$chisq,cancer_stratified_logrank_p=pchisq(slr$chisq,1,lower.tail=FALSE),
    unadjusted_HR=s0$conf.int['groupHigh','exp(coef)'],unadjusted_CI_low=s0$conf.int['groupHigh','lower .95'],
    unadjusted_CI_high=s0$conf.int['groupHigh','upper .95'],
    adjusted_HR=ss$conf.int['groupHigh','exp(coef)'],adjusted_CI_low=ss$conf.int['groupHigh','lower .95'],
    adjusted_CI_high=ss$conf.int['groupHigh','upper .95'],adjusted_p=ss$coefficients['groupHigh','Pr(>|z|)'],
    PH_group_p=ph['group','p'],PH_global_p=ph['GLOBAL','p'])
  check_equal(result,tab(paste0(s$panel,'/survival_comparison.tsv')))
  results[[s$panel]] <- result
  ph_rows[[s$panel]] <- data.frame(panel=s$panel,term=rownames(ph),ph,row.names=NULL)
  km <- survfit(Surv(months,OS_event)~group,data=d,conf.type='log-log')
  sk <- summary(km,censored=TRUE)
  curve <- data.frame(split=s$split,marker=s$marker,group=sub('group=','',as.character(sk$strata)),
    months=sk$time,survival=sk$surv,lower=sk$lower,upper=sk$upper,n_risk=sk$n.risk,n_event=sk$n.event,n_censor=sk$n.censor)
  rt <- summary(km,times=seq(0,120,24),extend=TRUE)
  risk <- data.frame(split=s$split,marker=s$marker,group=sub('group=','',as.character(rt$strata)),months=rt$time,n_risk=rt$n.risk)
  check_equal(curve,tab(paste0(s$panel,'/KM_curve_data.tsv')))
  curves[[s$panel]] <- curve; risks[[s$panel]] <- risk
  check_equal(risk,tab(paste0(s$panel,'/risk_table.tsv')))

}
stats <- do.call(rbind,results)
for(split in unique(specs$split)) {
  ix <- which(stats$split==split)
  for(name in c('logrank_p','cancer_stratified_logrank_p','adjusted_p')) stats[ix,paste0(name,'_BH')] <- p.adjust(stats[ix,name],'BH')
}
for(i in seq_len(nrow(specs))) check_equal(stats[i,,drop=FALSE],tab(paste0(specs$panel[i],'/survival_comparison.tsv')))



list(statistics = stats, curves = curves, risks = risks, PH = ph_rows, warnings = warnings)

}
