# Load unchanged CRAN DWLS 0.1.0 core functions without the unused Seurat
# signature-building dependency. No blocked DLL or security setting is touched.
load_dwls_core <- function(source_dir, optimized=FALSE) {
  env <- new.env(parent=baseenv())
  env$solve.QP <- quadprog::solve.QP
  env$lm <- stats::lm
  env$sd <- stats::sd
  for(name in c("solveOLSInternal","findDampeningConstant","solveDampenedWLSj","solveDampenedWLS"))
    sys.source(file.path(source_dir,paste0(name,".R")),envir=env)
  if(optimized) {
    fast_damping <- function(S,B,goldStandard) {
      ws <- as.vector((1/(S %*% goldStandard))^2)
      scaled <- ws/min(ws)
      finite_ws <- if(max(scaled)=="Inf") scaled[-which(scaled=="Inf")] else scaled
      subsets <- lapply(1:100,function(i){set.seed(i);sample(length(ws),size=length(ws)*.5)})
      solutions_sd <- NULL
      for(j in 1:ceiling(log2(max(finite_ws)))) {
        w <- scaled; w[which(w>2^(j-1))] <- 2^(j-1)
        solutions <- vapply(subsets,function(idx) {
          coef <- stats::lm.wfit(S[idx,,drop=FALSE],B[idx],w=w[idx])$coefficients
          coef*sum(goldStandard)/sum(coef)
        },numeric(ncol(S)))
        solutions_sd <- cbind(solutions_sd,apply(solutions,1,sd))
      }
      which.min(colMeans(solutions_sd^2))
    }
    fast_weighted <- function(S,B,goldStandard,j) {
      ws <- as.vector((1/(S %*% goldStandard))^2)
      w <- ws/min(ws); w[which(w>2^(j-1))] <- 2^(j-1)
      D <- crossprod(S,S*w); d <- crossprod(S,B*w)
      sc <- norm(D,"2")
      z <- quadprog::solve.QP(D/sc,d/sc,diag(ncol(S)),rep(0,ncol(S)))$solution
      names(z) <- colnames(S); z
    }
    environment(fast_damping) <- env; environment(fast_weighted) <- env
    env$findDampeningConstant <- fast_damping
    env$solveDampenedWLSj <- fast_weighted
  }
  env
}
