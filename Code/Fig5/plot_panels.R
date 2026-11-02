# Sourced by reproduce_Fig5.R after all numerical checks have passed.
source(file.path(code_dir, "shared/plotting.R"), local = TRUE, encoding = "UTF-8")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (panel in panels) {
  set.seed(123)
  warning_messages <- character()
  withCallingHandlers(
    if (panel %in% LETTERS[1:3]) {
      plot_files <- c(plot_files, plot_box_panel(panel))
    } else {
      source(file.path(code_dir, panel, "plot.R"), local = TRUE, encoding = "UTF-8")
    },
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warning_messages)) {
    message("Panel ", panel, ": ", length(warning_messages), " plotting warnings; original LOESS retained.")
    for (warning_text in unique(warning_messages)) message("  ", warning_text)
  }
}
cat("PDF/PNG output: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
