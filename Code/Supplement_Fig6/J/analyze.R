tab <- function(p) read_tsv(file.path(data_dir, p))
warnings <- character()
jd <- prepare_patients(tab('J/analysis_input.tsv'))
check_table(jd, context$patients, "patient", c("cancer","age","sex","stage","OS_time_days","OS_event","ImmuICBscore","CD45_TPM"))
groups <- unique(tab('shared/survival/lineage_mapping.tsv')$group)
ilr <- paste0('major_ilr',seq_len(length(groups)-1))
basis <- contr.helmert(length(groups)); basis <- sweep(basis,2,sqrt(colSums(basis^2)),'/')
rownames(basis) <- groups; colnames(basis) <- ilr
stopifnot(length(groups)==12,max(abs(crossprod(basis)-diag(11)))<1e-12)
jd[,ilr] <- log(as.matrix(jd[,groups])+1e-5) %*% basis
jd$score_z <- as.numeric(scale(jd$ImmuICBscore)); jd$cd45_z <- as.numeric(scale(log2(jd$CD45_TPM+1)))
base <- c('age10','sex','stage')
models <- list(Clinical_Score=c(base,'score_z'),Clinical_CD45_Score=c(base,'cd45_z','score_z'),
  Clinical_Major_Score=c(base,ilr,'score_z'),Clinical_Major_CD45_Score=c(base,ilr,'cd45_z','score_z'))
coefs <- list(); diagnostic <- list()
for(id in names(models)) {
  f <- as.formula(paste('Surv(OS_time_days,OS_event)~',paste(c(models[[id]],'strata(cancer)'),collapse='+')))
  fit <- withCallingHandlers(coxph(f,data=jd,ties='efron',x=TRUE,y=TRUE,model=TRUE,singular.ok=FALSE,
      control=coxph.control(iter.max=100)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  ss <- summary(fit)
  coefs[[id]] <- data.frame(model=id,term=rownames(ss$coefficients),beta=ss$coefficients[,'coef'],
    HR=ss$conf.int[,'exp(coef)'],CI_low=ss$conf.int[,'lower .95'],CI_high=ss$conf.int[,'upper .95'],
    p=ss$coefficients[,'Pr(>|z|)'],row.names=NULL)
  ph <- cox.zph(fit)$table
  diagnostic[[id]] <- data.frame(model=id,term=rownames(ph),ph,row.names=NULL)
}
all <- do.call(rbind,coefs); forest <- all[all$term=='score_z',]
reference <- tab('J/forest_plot_data.tsv')
check_equal(forest,reference)
max_error <- max(abs(as.matrix(forest[,3:7])-as.matrix(reference[,3:7])))

ref_all <- tab("J/reference_coefficients.tsv")
ref_all <- ref_all[ref_all$model %in% names(models), ]
check_table(all, ref_all, c("model", "term"), tolerance = 1e-8)
analysis$J <- list(data = forest, coefficients = all, PH = diagnostic, warnings = warnings)
