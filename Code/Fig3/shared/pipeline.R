# Stage orchestration shared by both main figures. Historical inputs are read-only.
pipeline_methods <- c("ImmuCellAI2", "BayesPrism", "CIBERSORT", "CITMIC",
                      "DWLS", "ImmuCellAI", "MuSiC_basic")
pipeline_labels <- setNames(c("ImmuCellAI 2.0", "BayesPrism", "CIBERSORT",
  "CITMIC", "DWLS", "ImmuCellAI", "MuSiC"), pipeline_methods)
pipeline_read <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM", na.strings = c("", "NA", "NaN", "nan"))
}
pipeline_matrix <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    m <- readRDS(path)
    if (is.list(m) && !is.null(m$prediction)) m <- m$prediction
  } else {
    first <- readLines(path, n = 1L, warn = FALSE)
    sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
    m <- as.matrix(read.table(path, header = TRUE, sep = sep, row.names = 1L,
      check.names = FALSE, comment.char = "", fileEncoding = "UTF-8-BOM"))
  }
  stopifnot(is.matrix(m), !is.null(rownames(m)), !is.null(colnames(m)),
            !anyDuplicated(rownames(m)), !anyDuplicated(colnames(m)))
  storage.mode(m) <- "double"
  m
}
pipeline_inputs <- function(panel) {
  panel <- match.arg(panel, if (figure_name == "Fig2") LETTERS[1:5] else LETTERS[1:5])
  cfg <- source(file.path(code_dir, panel, "input.R"), local = TRUE, encoding = "UTF-8")$value
  cfg$panel <- panel
  cfg$figure <- figure_name
  cfg
}
pipeline_statistics <- function(panel, run_dir = NULL, methods = pipeline_methods,
                                verify = is.null(run_dir)) {
  cfg <- pipeline_inputs(panel)
  if (cfg$kind %in% c("summary", "runtime")) {
    return(source(file.path(code_dir, panel, "statistics.R"), local = TRUE, encoding = "UTF-8")$value)
  }
  pipeline_correlations(cfg, run_dir, methods, verify)
}
pipeline_truth <- function(cfg) {
  if (cfg$kind == "legacy") {
    d <- pipeline_read(cfg$truth)
    return(data.frame(sample_id = d$sample_id, target = d$target_cell_type,
                      truth = d$target_fraction))
  }
  if (cfg$kind == "external") {
    m <- pipeline_matrix(cfg$truth)
    return(do.call(rbind, lapply(colnames(m), function(t) data.frame(
      sample_id = rownames(m), target = t, truth = m[, t]))))
  }
  # Frozen sample/target pairs supply truth only; predictions are reloaded below.
  d <- pipeline_read(cfg$pairs)
  d <- d[, c("sample_id", "target", "truth")]
  d <- unique(d)
  stopifnot(!anyDuplicated(d[c("sample_id", "target")]))
  if (cfg$kind == "npc") {
    denom <- pipeline_read(cfg$pairs)$denominator
    stopifnot(all(denom == "all_cells"))
  }
  d
}
pipeline_states <- function(cfg, method, target, columns) {
  if (cfg$kind %in% c("legacy", "external")) {
    env <- new.env(parent = baseenv())
    source(file.path(code_dir, cfg$panel, "mapping.R"), local = env, encoding = "UTF-8")
    name <- switch(method, CITMIC = "citmic_target_mapper",
      ImmuCellAI = "immucellai_target_mapper", CIBERSORT = "lm22_target_mapper",
      "state_target_mapper")
    states <- env[[name]](target)[[target]]
    # Historical runners evaluated the supported subset of an aggregate.
    return(intersect(states, columns))
  }
  map <- pipeline_read(cfg$mapping)
  row <- map[map$method == method & map$target == target, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  states <- if (is.na(row$predicted_states) || !nzchar(row$predicted_states)) character() else
    trimws(strsplit(row$predicted_states, "|", fixed = TRUE)[[1]])
  stopifnot(all(states %in% columns))
  states
}
pipeline_prediction_files <- function(cfg, run_dir) {
  if (is.null(run_dir)) return(cfg$predictions)
  root <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  lock_path <- file.path(root, "pipeline_job.rds")
  if (!file.exists(lock_path)) stop("Not a pipeline run: missing pipeline_job.rds")
  lock <- readRDS(lock_path)
  stopifnot(identical(lock$figure, cfg$figure), identical(lock$panel, cfg$panel),
    identical(lock$source_hashes, unname(tools::md5sum(lock$source_files))))
  relative <- cfg$prediction_relative
  setNames(file.path(root, relative), names(relative))
}
pipeline_correlations <- function(cfg, run_dir, methods, verify) {
  stopifnot(length(methods) > 0L, !anyDuplicated(methods), all(methods %in% pipeline_methods))
  files <- pipeline_prediction_files(cfg, run_dir)
  missing <- methods[!file.exists(files[methods])]
  if (length(missing)) stop("Missing sample-level predictions for ", cfg$figure, " ", cfg$panel,
    ": ", paste(missing, collapse = ", "),
    ". Historical summary values are NOT substituted. Use plot() for the archived figure.")
  truth <- pipeline_truth(cfg)
  expected <- if (figure_name == "Fig2") pipeline_read(cfg$plot_data) else
    read.csv(cfg$plot_data, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  targets <- if (figure_name == "Fig2") unique(expected$target) else unique(expected$cellType)
  if (!is.null(cfg$target_alias)) targets <- unname(cfg$target_alias[targets])
  stopifnot(!anyNA(targets), all(targets %in% truth$target))
  ids <- unique(truth$sample_id)
  rows <- pairs <- list()
  for (method in methods) {
    pred <- pipeline_matrix(files[[method]])
    if (!setequal(ids, rownames(pred)) && setequal(ids, colnames(pred))) pred <- t(pred)
    if (!setequal(ids, rownames(pred))) stop("Sample IDs differ from frozen truth: ", method)
    pred <- pred[ids, , drop = FALSE]
    if (any(!is.finite(pred))) stop("Non-finite raw predictions: ", method)
    for (target in targets) {
      d <- truth[truth$target == target, , drop = FALSE]
      states <- pipeline_states(cfg, method, target, colnames(pred))
      y <- if (length(states)) rowSums(pred[d$sample_id, states, drop = FALSE]) else
        rep(NA_real_, nrow(d))
      ok <- is.finite(d$truth) & is.finite(y)
      status <- if (!length(states)) "unavailable" else if (sum(ok) < 3L)
        "undefined_correlation" else if (stats::sd(y[ok]) == 0)
        "constant_prediction" else if (stats::sd(d$truth[ok]) == 0)
        "undefined_correlation" else "valid"
      r <- if (status == "valid") stats::cor(d$truth[ok], y[ok]) else NA_real_
      label <- if (is.null(cfg$target_alias)) target else names(cfg$target_alias)[match(target, cfg$target_alias)]
      rows[[length(rows) + 1L]] <- data.frame(method = method, target = label,
        n_designed = nrow(d), n_pairs = sum(ok), Pearson = r, status = status)
      pairs[[length(pairs) + 1L]] <- data.frame(method = method, target = label,
        sample_id = d$sample_id, truth = d$truth, prediction = y)
    }
  }
  correlations <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(methods, function(m) {
    r <- correlations$Pearson[correlations$method == m]
    r <- r[is.finite(r)]
    data.frame(method = m, n_valid = length(r), n_targets = length(targets),
      mean = if (length(r)) mean(r) else NA_real_,
      SD = if (length(r) > 1L) stats::sd(r) else NA_real_)
  }))
  if (verify) {
    if (figure_name == "Fig3") {
      expected$method <- names(pipeline_labels)[match(expected$method, pipeline_labels)]
      expected$target <- expected$cellType
      expected$Pearson <- as.numeric(expected$correlation)
    }
    expected <- expected[expected$method %in% methods, ]
    key <- function(x) paste(x$method, x$target, sep = "\r")
    stopifnot(!anyDuplicated(key(correlations)), !anyDuplicated(key(expected)),
              setequal(key(correlations), key(expected)))
    ref <- expected$Pearson[match(key(correlations), key(expected))]
    stopifnot(identical(is.na(ref), is.na(correlations$Pearson)),
      all(abs(ref[is.finite(ref)] - correlations$Pearson[is.finite(ref)]) < 1e-10))
  }
  list(figure = figure_name, panel = cfg$panel, correlations = correlations,
       summary = summary, pairs = do.call(rbind, pairs), verified_against_archive = verify)
}
pipeline_predict <- function(panel, confirm = FALSE, cores = 2L) {
  cfg <- pipeline_inputs(panel)
  if (cfg$kind %in% c("summary", "runtime"))
    stop("This panel does not have a new cohort-level seven-tool prediction stage.")
  if (!isTRUE(confirm)) stop("Prediction is expensive. Call predict(panel, confirm = TRUE) explicitly.")
  stopifnot(length(cores) == 1L, is.finite(cores), cores >= 1, cores == as.integer(cores))
  pipeline_run(cfg, cores)
}
pipeline_run <- function(cfg, cores) {
  source_files <- unique(c(cfg$bulk, cfg$truth, cfg$pairs, cfg$mapping, cfg$reference, cfg$markers,
                           cfg$design, cfg$locked_protocol, cfg$baseline))
  source_files <- source_files[!is.na(source_files) & nzchar(source_files)]
  stopifnot(all(file.exists(source_files)))
  parent <- file.path(code_dir, "output_analysis")
  dir.create(parent, showWarnings = FALSE)
  root <- tempfile(paste0(cfg$panel, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = parent)
  stopifnot(!file.exists(root))
  dir.create(root)
  job <- cfg
  job$output <- normalizePath(root, winslash = "/", mustWork = TRUE)
  job$source_files <- source_files
  job$source_hashes <- unname(tools::md5sum(source_files))
  job$code_dir <- code_dir
  job$cores <- as.integer(cores)
  job$library_root <- getOption("Immu2.library_root",
    "C:/Users/17729/Documents/Codex/2026-05-26/files-mentioned-by-the-user-immunehierdeconv")
  job$engine <- file.path(code_dir, cfg$engine)
  job_file <- file.path(root, "pipeline_job.rds")
  saveRDS(job, job_file)
  steps <- if (cfg$kind %in% c("bg", "npc")) c("prepare", pipeline_methods) else
    if (cfg$kind %in% c("legacy", "external")) "all" else pipeline_methods
  runner <- file.path(code_dir, "shared/prediction_worker.R")
  for (task in steps) {
    message(cfg$figure, " ", cfg$panel, ": ", task)
    rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    status <- system2(rscript,
      c(shQuote(runner), shQuote(job_file), shQuote(task)))
    if (!identical(status, 0L)) stop("Prediction stopped at ", task, ". Partial results: ", root)
  }
  result <- pipeline_statistics(cfg$panel, run_dir = root, verify = FALSE)
  # Necessary analysis result, not an automatic checksum/validation report.
  saveRDS(result, file.path(root, "statistics.rds"))
  message("Predictions and statistics: ", root)
  invisible(list(run_dir = root, statistics = result))
}
pipeline_plot <- function(results = NULL, output_dir = NULL,
                          format = "both", dpi = 300L) {
  if (is.null(output_dir)) output_dir <- if (is.null(results)) file.path(code_dir, "output_plot") else
    file.path(code_dir, "output_plot", paste0("recomputed_", base::format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (is.null(results)) {
    old <- options(setNames(list("--plots", output_dir, format, dpi),
      paste0(figure_name, c(".args", ".output_dir", ".format", ".dpi"))))
    on.exit(options(old), add = TRUE)
    return(source(file.path(code_dir, paste0("reproduce_", figure_name, ".R")),
      local = new.env(parent = globalenv()), encoding = "UTF-8")$value)
  }
  # New results must include every cohort. Never silently mix an absent cohort with cached r.
  required <- if (figure_name == "Fig2") LETTERS[1:4] else LETTERS[1:3]
  stopifnot(is.list(results), setequal(names(results), required))
  analysis <- list()
  data_dir <- file.path(project, "Data", figure_name)
  source(file.path(code_dir, "shared/statistics.R"), local = TRUE, encoding = "UTF-8")
  for (panel in required) {
    cfg <- pipeline_inputs(panel)
    x <- results[[panel]]
    stopifnot(x$figure == figure_name, x$panel == panel)
    d <- if (figure_name == "Fig2") pipeline_read(cfg$plot_data) else
      read.csv(cfg$plot_data, check.names = FALSE, fileEncoding = "UTF-8-BOM")
    key <- function(m, t) paste(m, t, sep = "\r")
    k <- if (figure_name == "Fig2") key(d$method, d$target) else
      key(names(pipeline_labels)[match(d$method, pipeline_labels)], d$cellType)
    actual <- key(x$correlations$method, x$correlations$target)
    stopifnot(!anyDuplicated(actual), setequal(actual, k))
    c <- x$correlations[match(k, actual), ]
    d$status <- c$status
    if (figure_name == "Fig2") {
      d$Pearson <- c$Pearson
      d$sector_fraction <- abs(c$Pearson)
      d$symbol <- ifelse(c$status == "valid", "sector", ifelse(c$status == "unavailable", "/", "x"))
      d$n <- c$n_designed
    } else {
      d$r <- d$correlation <- c$Pearson
    }
    analysis[[panel]] <- d
  }
  rows <- list()
  for (panel in required) for (m in pipeline_methods) {
    d <- analysis[[panel]]
    keep <- if (figure_name == "Fig2") d$method == m else d$method == pipeline_labels[[m]]
    r <- if (figure_name == "Fig2") d$Pearson[keep] else d$r[keep]
    valid <- is.finite(r); n <- sum(valid)
    if (n < 2L) stop("At least two evaluable targets required for the existing SD renderer: ", panel, "/", m)
    avg <- mean(r[valid]); std <- stats::sd(r[valid])
    rows[[length(rows) + 1L]] <- if (figure_name == "Fig2") data.frame(
      dataset = d$dataset[1], scope = "all_valid", method = m, n = n, n_candidate = length(r),
      mean_Pearson = avg, SD_Pearson = std, lower = avg-std, upper = avg+std) else data.frame(
      Dataset = d$Dataset[1], source_panel = panel, method = pipeline_labels[[m]], n_valid = n,
      n_candidates = length(r), mean = avg, SD = std, lower = avg-std, upper = avg+std,
      included_targets = paste(d$cellType[keep][valid], collapse = " | "))
  }
  analysis[[if (figure_name == "Fig2") "E" else "D"]] <- do.call(rbind, rows)
  if (figure_name == "Fig3") {
    # Runtime is the separately measured historical experiment, not retimed predictions.
    e <- new.env(parent = environment())
    source(file.path(code_dir, "shared/statistics.R"), local = e, encoding = "UTF-8")
    e$analysis <- analysis
    source(file.path(code_dir, "E/analyze.R"), local = e, encoding = "UTF-8")
    analysis$E <- e$analysis$E
  }
  out <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  root <- tolower(paste0(normalizePath(code_dir, winslash = "/"), "/"))
  if (!startsWith(tolower(paste0(out, "/")), root))
    stop("New-result figures must be saved under this Code figure directory.")
  stopifnot(format %in% c("pdf", "png", "both"), dpi >= 72L)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  output_dir <- out; output_format <- format; panels <- LETTERS[1:5]; plot_files <- character()
  source(file.path(code_dir, "plot_panels.R"), local = TRUE, encoding = "UTF-8")
  invisible(list(analysis = analysis, files = plot_files))
}
