# Summary from recomputed panel statistics; no silent fallback to archived correlations.
if (!is.null(run_dir)) stop("For new runs, pass named A-D/A-C result lists to plot(results).")
do.call(rbind, lapply(cfg$panels, function(panel) {
  x <- pipeline_statistics(panel, methods = methods, verify = verify)$summary
  cbind(panel = panel, x)
}))
