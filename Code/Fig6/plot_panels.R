# Sourced by reproduce_Fig6.R after all numerical checks have passed.
source(file.path(code_dir, "shared/plotting.R"), local = TRUE, encoding = "UTF-8")
for (panel in panels) {
  warning_messages <- character()
  withCallingHandlers(
    if (panel %in% c("B", "C", "D")) {
      plot_files <- c(plot_files, plot_comparison_panel(panel))
    } else {
      source(file.path(code_dir, panel, "plot.R"), local = TRUE, encoding = "UTF-8")
    },
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warning_messages)) {
    message("Panel ", panel, ": ", length(warning_messages), " plotting warnings.")
    for (warning_text in unique(warning_messages)) message("  ", warning_text)
  }
}
cat("PDF/PNG output: ", output_dir, "\n", sep = "")
