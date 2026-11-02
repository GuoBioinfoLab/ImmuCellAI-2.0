# Optional fixed-protocol reproduction, not cohort/parameter/seed selection.
local({
source_files <- lapply(sys.frames(), function(frame) frame$ofile)
source_files <- Filter(function(path) is.character(path) && length(path) == 1L &&
                         !is.na(path) && nzchar(path), source_files)
if (length(source_files)) {
  script_file <- source_files[[length(source_files)]]
} else {
  file_arg <- grep("^--file=", commandArgs(), value=TRUE)
  if (!length(file_arg)) stop("Run this file with source() or Rscript.")
  script_file <- sub("^--file=", "", file_arg[1L])
}
code <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
input <- normalizePath(file.path(code, "../../Data/Fig4"), winslash="/")
suppressPackageStartupMessages(library(ranger))
read_tab <- function(p) read.delim(file.path(input,p), check.names=FALSE, stringsAsFactors=FALSE)
arch <- read_tab("H/predictions.tsv")
arch <- arch[order(arch$sample), ]
fine <- read_tab("H/fine_fractions.tsv")
canonical <- read_tab("H/canonical_fractions.tsv")
saved <- read_tab("H/ROC_scores.tsv")
stopifnot(nrow(arch)==334, !anyDuplicated(arch$sample), setequal(arch$fold,1:5))
align <- function(d) {
  stopifnot(!anyDuplicated(d$sample), setequal(d$sample,arch$sample))
  d <- d[match(arch$sample,d$sample), ]
  stopifnot(identical(d$sample,arch$sample))
  d
}
fine <- align(fine); canonical <- align(canonical); saved <- align(saved)
stopifnot(identical(arch$ResponseBinary,saved$ResponseBinary))
clr <- function(d) {
  x <- as.matrix(d[,setdiff(names(d),"sample"),drop=FALSE])
  stopifnot(all(is.finite(x)), all(x>=0), max(abs(rowSums(x)-1))<1e-8)
  lx <- log(pmax(x,0)+1e-5)
  sweep(lx,1,rowMeans(lx),"-")
}
features <- list(Fine53_RF=clr(fine),Canonical12_RF=clr(canonical))
stopifnot(ncol(features$Fine53_RF)==53,ncol(features$Canonical12_RF)==12)
prediction <- data.frame(sample=arch$sample,fold=arch$fold,ResponseBinary=arch$ResponseBinary)
errors <- list()
for (model in names(features)) {
  message("Reproducing fixed five-fold model: ",model)
  x <- features[[model]]
  pp <- rep(NA_real_,nrow(x))
  for (fold in 1:5) {
    tr <- arch$fold!=fold; va <- !tr
    center <- colMeans(x[tr,,drop=FALSE])
    scalev <- apply(x[tr,,drop=FALSE],2,sd)
    scalev[!is.finite(scalev)|scalev==0] <- 1
    transform <- function(z) as.data.frame(sweep(sweep(z,2,center,"-"),2,scalev,"/"),check.names=FALSE)
    yy <- arch$ResponseBinary[tr]
    seed <- 9001L+fold
    set.seed(seed)
    fit <- ranger(x=transform(x[tr,,drop=FALSE]),y=factor(yy,levels=c(0,1)),
                  probability=TRUE,classification=TRUE,num.trees=1200,
                  mtry=max(1,floor(sqrt(ncol(x)))),min.node.size=8,
                  class.weights=c("0"=0.5/mean(yy==0),"1"=0.5/mean(yy==1)),
                  seed=seed,num.threads=8)
    pp[va] <- predict(fit,data=transform(x[va,,drop=FALSE]))$predictions[,"1"]
  }
  prediction[[model]] <- pp
  errors[[model]] <- max(abs(pp-saved[[model]]))
  stopifnot(all(is.finite(pp)),errors[[model]]<1e-10)
}
print(unlist(errors))
cat("Both RF models reproduce archived predictions within 1e-10.\n")
})
