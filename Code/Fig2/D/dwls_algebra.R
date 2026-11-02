# Algebraically equivalent speed-ups for native DWLS 0.1.0 internals.
# Keep the public solver, 100 bootstrap subsets, damping search and QP unchanged.
patch_dwls_algebra <- function() {
  ns <- asNamespace("DWLS")
  fast_damping <- function(S, B, goldStandard) {
    ws <- as.vector((1/(S %*% goldStandard))^2)
    wsScaled <- ws/min(ws)
    finite_ws <- if (max(wsScaled) == "Inf") wsScaled[-which(wsScaled == "Inf")] else wsScaled
    subsets <- lapply(1:100, function(i) {
      set.seed(i)
      sample(length(ws), size = length(ws)*0.5)
    })
    solutionsSd <- NULL
    for (j in 1:ceiling(log2(max(finite_ws)))) {
      w <- wsScaled
      w[which(wsScaled > 2^(j-1))] <- 2^(j-1)
      solutions <- vapply(subsets, function(idx) {
        coef <- stats::lm.wfit(S[idx, , drop=FALSE], B[idx], w=w[idx])$coefficients
        coef*sum(goldStandard)/sum(coef)
      }, numeric(ncol(S)))
      solutionsSd <- cbind(solutionsSd, apply(solutions,1,sd))
    }
    which.min(colMeans(solutionsSd^2))
  }
  fast_weighted <- function(S, B, goldStandard, j) {
    ws <- as.vector((1/(S %*% goldStandard))^2)
    w <- ws/min(ws)
    w[which(w > 2^(j-1))] <- 2^(j-1)
    D <- crossprod(S, S*w)
    d <- crossprod(S, B*w)
    sc <- norm(D, "2")
    result <- quadprog::solve.QP(D/sc, d/sc, diag(ncol(S)), rep(0,ncol(S)))$solution
    names(result) <- colnames(S)
    result
  }
  for (entry in list(list("findDampeningConstant",fast_damping),list("solveDampenedWLSj",fast_weighted))) {
    environment(entry[[2]]) <- ns
    locked <- bindingIsLocked(entry[[1]],ns)
    if(locked) unlockBinding(entry[[1]],ns)
    assign(entry[[1]],entry[[2]],envir=ns)
    if(locked) lockBinding(entry[[1]],ns)
  }
  invisible(TRUE)
}
