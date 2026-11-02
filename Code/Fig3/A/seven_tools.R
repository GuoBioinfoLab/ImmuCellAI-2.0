# Execution is controlled by shared/prediction_worker.R; source() alone is not a run.
if (!exists("job", inherits = FALSE)) stop("Use the run_pipeline.R API to run predictions.")
# Fixed seven-method comparison on measured BostonGene blood data.
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
options(stringsAsFactors = FALSE, warn = 1, digits = 17)
workspace <- job$library_root; cores <- job$cores
.libPaths(c(file.path(workspace, "Rlib_GSE176078_comparison"), file.path(workspace, "Rlib_tools"),
            file.path(workspace, "Rlib_test"), .libPaths()))
Sys.setenv(LC_ALL = "English_United States.utf8", LANG = "English_United States.utf8",
           R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
script <- job$engine
out <- job$output
cohort <- job$cohort; batch <- dirname(cohort)
response <- dirname(dirname(batch))
helper <- file.path(job$code_dir, "A/legacy_helpers.R")
dwls.loader <- file.path(job$code_dir, "A/load_dwls_core.R")
dwls.source <- file.path(job$code_dir, "A/dwls_core")
source(helper); run_immucellai <- TRUE
methods <- c("ImmuCellAI2", "BayesPrism", "CIBERSORT", "CITMIC", "DWLS", "ImmuCellAI", "MuSiC_basic")
stopifnot(task %in% c("prepare", methods), cores >= 1L)
for (d in c("input", "references", "predictions", "metrics", "logs"))
  dir.create(file.path(out, d), showWarnings = FALSE, recursive = TRUE)
md5 <- function(p) unname(tools::md5sum(p))
write_json <- function(x, path) jsonlite::write_json(x, file.path(out, path), pretty = TRUE,
                                                    auto_unbox = TRUE, digits = 16, na = "null")
write_tsv <- function(x, path) fwrite(x, file.path(out, path), sep = "\t", na = "NA")
write_mat <- function(x, path, id = "sample_id") {
  d <- data.table(ID = rownames(x), x, check.names = FALSE); setnames(d, "ID", id)
  write_tsv(d, path)
}
read_mat <- function(path) {
  d <- fread(path, check.names = FALSE)
  x <- as.matrix(d[, -1, with = FALSE]); rownames(x) <- d[[1]]; storage.mode(x) <- "double"
  stopifnot(!anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)))
  x
}
original.protocol <- jsonlite::read_json(file.path(batch, "locked_protocol.json"), simplifyVector = FALSE)
original.entries <- original.protocol$cohorts$BG_blood$evaluation
entries <- Filter(function(e) length(e$states) > 0L, original.entries)
targets <- vapply(entries, `[[`, character(1), "target")
excluded <- c("Th1_17", "cNK", "NKreg", "Eosinophil", "Basophil", "Neutrophil", "CD4T", "CD8T")
selected <- setdiff(targets, excluded)
stopifnot(length(targets) == 34L, length(selected) == 26L)
source.files <- c(file.path(cohort, "bulk_linear.tsv.gz"), file.path(cohort, "truth_published.tsv"),
                  file.path(cohort, "results/ImmuCellAI2_53states.tsv"),
                  file.path(cohort, "results/per_cell_pearson.tsv"), file.path(cohort, "results/settings.json"),
                  file.path(batch, "locked_protocol.json"))
source.hash <- as.list(setNames(md5(source.files), basename(source.files)))
protocol.path <- file.path(out, "input/locked_protocol.json")
cache <- file.path(out, "input/prepared.rds")

if (task == "prepare") {
  stopifnot(!file.exists(protocol.path))
  stopifnot(jsonlite::read_json(file.path(cohort, "results/independent_verification.json"))$passed)
  bulk <- read_mat(source.files[1]); truth <- fread(source.files[2], check.names = FALSE)
  baseline <- read_mat(source.files[3])
  settings <- jsonlite::read_json(source.files[5], simplifyVector = TRUE)
  stopifnot(ncol(bulk) == 45L, identical(colnames(bulk), truth$sample_id),
            identical(rownames(baseline), truth$sample_id), all(is.finite(bulk)), all(bulk >= 0),
            settings$n.iter == 50, settings$hierarchy.mode == "tcell", !settings$add.unknown)
  old.hash <- jsonlite::read_json(file.path(cohort, "results/input_md5.json"), simplifyVector = TRUE)
  stopifnot(old.hash$bulk == md5(source.files[1]), old.hash$truth == md5(source.files[2]),
            old.hash$protocol == md5(source.files[6]))
  ref.full <- as.matrix(ImmuCellAI2.0::load_immucellai2_reference())
  provenance <- jsonlite::read_json(file.path(batch, "model_provenance.json"), simplifyVector = TRUE)
  stopifnot(md5(provenance$reference_file) == provenance$reference_md5)
  genes <- intersect(rownames(ref.full), rownames(bulk))
  genes <- genes[rowSums(ref.full[genes, , drop = FALSE]) > 0]
  ref <- ref.full[genes, , drop = FALSE]
  signature <- select_signature_genes(ref, 1000)
  e <- new.env(); data("LM22", package = "CIBERSORT", envir = e); lm22 <- as.matrix(e$LM22)
  expected <- list(ImmuCellAI2 = colnames(ref), BayesPrism = colnames(ref), DWLS = colnames(ref), MuSiC_basic = colnames(ref),
                   CIBERSORT = colnames(lm22), CITMIC = colnames(CITMIC::GetData_CITMIC("matrix_cell_go_inter")),
                   ImmuCellAI = c("DC", "Bcell", "Monocyte", "Macrophage", "NK", "Neutrophil", "CD4_T", "CD8_T", "NKT", "Gamma_delta",
                                  "CD4_naive", "Tr1", "nTreg", "iTreg", "Th1", "Th2", "Th17", "Tfh", "CD8_naive", "Cytotoxic", "Exhausted", "MAIT", "Central_memory", "Effector_memory"))
  native <- list(
    CIBERSORT = list(B = c("B cells naive", "B cells memory", "Plasma cells"), Bnaive = "B cells naive",
      Nonplasma_B = c("B cells naive", "B cells memory"), Plasma_lineage = "Plasma cells",
      CD4T = c("T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated", "T cells follicular helper", "T cells regulatory (Tregs)"),
      CD8T = "T cells CD8", CD4Tnaive = "T cells CD4 naive", Tfh = "T cells follicular helper", Treg = "T cells regulatory (Tregs)",
      NK = c("NK cells resting", "NK cells activated"), Monocyte = "Monocytes", Eosinophil = "Eosinophils", Neutrophil = "Neutrophils",
      DC_total = c("Dendritic cells resting", "Dendritic cells activated")),
    CITMIC = list(B = "B cells", Bnaive = "Naive B cell", Plasma_lineage = "Plasma cells",
      CD4T = "CD4+ T cells", CD8T = "CD8+ T cells", CD4Tnaive = "Naive CD4+ T cell", CD8Tnaive = "Naive CD8+ T cell",
      CD4Tcm = "CD4+ Tcm", CD4Tem = "CD4+ Tem", CD8Tcm = "CD8+ Tcm", CD8Tem = "CD8+ Tem",
      Th1 = "Th1 cells", Th2 = "Th2 cells", Tfh = "Tfh", Treg = "Tregs", NK = "NK cells", Monocyte = "Monocytes",
      Neutrophil = "Neutrophils", Eosinophil = "Eosinophils", Basophil = "Basophils", pDC = "pDCs", cDC = "cDCs", DC_total = "DCs", T_total = "T cells"),
    ImmuCellAI = list(B = "Bcell", CD4T = "CD4_T", CD8T = "CD8_T", CD4Tnaive = "CD4_naive", CD8Tnaive = "CD8_naive",
      Th1 = "Th1", Th2 = "Th2", Tfh = "Tfh", Treg = c("nTreg", "iTreg"), NK = "NK", Monocyte = "Monocyte", Neutrophil = "Neutrophil", DC_total = "DC"))
  mapping <- rbindlist(lapply(methods, function(method) rbindlist(lapply(entries, function(entry) {
    states <- if (method %in% c("ImmuCellAI2", "BayesPrism", "DWLS", "MuSiC_basic")) unlist(entry$states) else native[[method]][[entry$target]]
    if (is.null(states)) states <- character()
    stopifnot(!anyDuplicated(states), all(states %in% expected[[method]]))
    note <- if (!length(states)) "No sufficiently specific native output; no broader lineage substitution or overlapping score subtraction."
      else if (method %in% c("ImmuCellAI2", "BayesPrism", "DWLS", "MuSiC_basic")) entry$note
      else "Native label-based mapping; flow gates and signature definitions may differ. Aggregates are explicit, not interchangeable with subtypes."
    if (method %in% c("CIBERSORT", "CITMIC") && entry$target == "Plasma_lineage") note <- "Native plasma-cell signature as plasma-lineage correspondence; does not separately identify plasmablasts."
    data.table(method = method, target = entry$target, group = entry$group, primary = entry$primary,
                selected_display = entry$target %in% selected, predicted_states = paste(states, collapse = " | "),
                truth_terms = paste(unlist(entry$truth_terms), collapse = " | "),
                mapping_kind = if (length(states)) "operational_label_match" else "unavailable", mapping_note = note)
  }))))
  packages <- c(ImmuCellAI2 = "ImmuCellAI2.0", BayesPrism = "BayesPrism", CIBERSORT = "CIBERSORT", CITMIC = "CITMIC", DWLS = "DWLS", ImmuCellAI = "ImmuCellAI", MuSiC_basic = "MuSiC")
  parameters <- list(ImmuCellAI2 = settings,
    BayesPrism = list(chain.length = 600, burn.in = 500, thinning = 2, alpha = 1, seed = 123, update.gibbs = FALSE, theta = "first/state"),
    CIBERSORT = list(signature = "native LM22", QN = FALSE, perm = 0, maxSize = 500),
    CITMIC = list(input = "log2(published TPM+1)", weighted = TRUE, base = 10, damping = .9),
    DWLS = list(reference_only_genes = 1000, solver = "solveDampenedWLS", optimization = "validated algebraic speedup"),
    ImmuCellAI = list(data_type = "RNA-seq", group_tag = FALSE, response_tag = FALSE, customer = 0),
    MuSiC_basic = list(S = 1, Sigma = 1e-8, iter.max = 1000, nu = 1e-4, eps = .01))
  write_tsv(mapping, "input/locked_mapping.tsv")
  write_mat(lm22, "references/CIBERSORT_LM22.tsv", "gene")
  writeLines(signature, file.path(out, "references/DWLS_genes.txt"))
  saveRDS(list(bulk = bulk, ref = ref, lm22 = lm22, signature = signature), cache)
  protocol <- list(dataset = "BostonGene blood", n_samples = 45, sample_ids = colnames(bulk),
    targets = targets, selected_targets = selected, excluded_from_display = excluded,
    source_files = as.list(source.files), source_md5 = source.hash, parameters = parameters,
    package_versions = lapply(packages, function(p) as.character(packageVersion(p))),
    reference_source = provenance$reference_file, reference_md5 = provenance$reference_md5,
    LM22_md5 = md5(file.path(out, "references/CIBERSORT_LM22.tsv")),
    mapping_md5 = md5(file.path(out, "input/locked_mapping.tsv")), cache_md5 = md5(cache),
    helper_paths = c(helper, dwls.loader), helper_md5 = as.list(setNames(md5(c(helper, dwls.loader)), c("legacy", "dwls"))),
    reference_policy = "Native LM22, CITMIC and ImmuCellAI signatures; original fixed 53-state atlas for the other four tools.",
    input_policy = "Same original published linear TPM matrix; native per-method transformations only. No sample selection and no marker tuning against flow.",
    metrics = "Paired finite observations per target; Pearson against unchanged published flow percentages / 100. Unavailable and constant predictions remain NA.",
    caveats = c("CIBERSORT R + LM22, not CIBERSORTx.", "MuSiC basic with static reference and constant variance, not full multi-donor MuSiC.",
      "BayesPrism receives continuous TPM as pseudo-counts and a static reference; settings are the existing comparison settings, not an asserted default raw-count workflow.",
      "CITMIC and ImmuCellAI return scores, not calibrated cell fractions. Do not renormalize or compare their MAE as fraction errors.",
      "2.0 n.iter=50 is retained after examining this cohort's iteration sensitivity; no claim of converged or independently tuned superiority.",
      "The 26-target display was selected after inspecting 2.0 results. Full 34-target and original primary10 results are retained for every method.",
      "Different tools cover different targets; own-coverage means are descriptive and not an unbiased ranking. Parent and subtype targets overlap.",
      "Mixed blood preparations and per-panel missing flow data remain unchanged. Original dataset limitations and reference-overlap uncertainty apply."),
    prepared_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  write_json(protocol, "input/locked_protocol.json")
  cat("PREPARED: 45 real samples, 34 full targets, 26 displayed targets, 238 locked mappings\n")
  quit(status = 0)
}

protocol <- jsonlite::read_json(protocol.path, simplifyVector = FALSE)
stopifnot(identical(protocol$source_md5, source.hash), md5(cache) == protocol$cache_md5,
          md5(file.path(out, "input/locked_mapping.tsv")) == protocol$mapping_md5,
          md5(protocol$reference_source) == protocol$reference_md5)
protocol.hash <- md5(protocol.path)
mapping <- fread(file.path(out, "input/locked_mapping.tsv"), na.strings = "NA")
parts <- function(x) if (is.na(x) || !nzchar(x)) character() else strsplit(x, " | ", fixed = TRUE)[[1]]

final <- file.path(out, "predictions", paste0(task, ".rds"))
if (file.exists(final)) {
  stopifnot(readRDS(final)$protocol_md5 == protocol.hash)
  cat("Verified completed method", task, "\n"); quit(status = 0)
}
d <- readRDS(cache); bulk <- d$bulk; ref <- d$ref; b <- bulk[rownames(ref), , drop = FALSE]
set.seed(123); started <- Sys.time(); details <- list()
if (task == "ImmuCellAI2") {
  pred <- read_mat(source.files[3])
  kind <- "immune_relative_fraction"
  details <- list(reused_verified_n_iter_50 = TRUE, source = source.files[3], source_md5 = md5(source.files[3]))
} else if (task == "BayesPrism") {
  suppressPackageStartupMessages(library(BayesPrism))
  h <- ImmuCellAI2.0::create_default_53_hierarchy(colnames(ref), hierarchy.mode = "tcell")
  h <- h[match(colnames(ref), h$state_name), ]
  prism <- BayesPrism::new.prism(reference = t(ref), input.type = "count.matrix", cell.type.labels = h$major_lineage,
    cell.state.labels = colnames(ref), key = NA_character_, mixture = t(b), outlier.cut = .01, outlier.fraction = .1, pseudo.min = 1e-8)
  fit <- BayesPrism::run.prism(prism, n.cores = cores, update.gibbs = FALSE,
    gibbs.control = list(chain.length = 600, burn.in = 500, thinning = 2, alpha = 1, seed = 123))
  pred <- BayesPrism::get.fraction(fit, which.theta = "first", state.or.type = "state")
  details <- list(genes_used = colnames(prism@mixture)); kind <- "static_reference_TPM_pseudocount_fraction"
} else if (task == "CITMIC") {
  pred <- CITMIC::CITMIC(as.data.frame(log2(bulk + 1)), weighted = TRUE, base = 10, damping = .9, cl.cores = cores)
  kind <- "native_enrichment_score_not_fraction"
} else if (task == "ImmuCellAI") {
  suppressPackageStartupMessages(library(ImmuCellAI)); patch_immucellai_gsva()
  fit <- ImmuCellAI::ImmuCellAI_new(sample = as.data.frame(bulk), data_type = "RNA-seq", group_tag = FALSE, response_tag = FALSE, customer = 0)
  pred <- fit$Sample_abundance; kind <- "native_abundance_score_not_fraction"
} else {
  if (task == "DWLS") {
    source(dwls.loader)
    original <- load_dwls_core(dwls.source, FALSE); fast <- load_dwls_core(dwls.source, TRUE)
    S <- ref[d$signature, , drop = FALSE]; B <- b[d$signature, , drop = FALSE]
    common.scale <- max(S); S <- S / common.scale; B <- B / common.scale
    tests <- c(1L, ncol(B))
    errors <- vapply(tests, function(j) {
      set.seed(123 + j); a <- original$solveDampenedWLS(S, B[, j])
      set.seed(123 + j); z <- fast$solveDampenedWLS(S, B[, j])
      max(abs(a - z))
    }, numeric(1))
    stopifnot(max(errors) < 1e-8)
    write_json(list(passed = TRUE, sample_indices = tests, max_abs_errors = errors), "logs/DWLS_algebra_check.json")
    worker <- function(j) { set.seed(123 + j); as.numeric(fast$solveDampenedWLS(S, B[, j])) }
    export <- c("S", "B", "fast"); columns <- colnames(S)
    details <- list(genes = d$signature, algebra_error = errors, common_scale = common.scale)
    kind <- "immune_reference_relative_fraction"
  } else if (task == "CIBERSORT") {
    genes <- intersect(rownames(d$lm22), rownames(bulk))
    S <- d$lm22[genes, , drop = FALSE]; B <- bulk[genes, , drop = FALSE]
    worker <- function(j) {
      z <- B[, c(j, j), drop = FALSE]; colnames(z) <- c("sample", "duplicate")
      stopifnot(max(z) >= 50)
      fit <- CIBERSORT::cibersort(S, z, perm = 0, QN = FALSE, maxSize = 500)
      as.numeric(fit[1, colnames(S)])
    }
    export <- c("S", "B"); columns <- colnames(S)
    details <- list(genes = genes, signature = "native LM22", implementation = "CIBERSORT R, not CIBERSORTx")
    kind <- "LM22_immune_relative_fraction"
  } else if (task == "MuSiC_basic") {
    X <- sweep(ref, 2, colSums(ref), "/"); Y <- sweep(b, 2, colSums(b), "/")
    S <- rep(1, ncol(X)); Sigma <- matrix(1e-8, nrow(X), ncol(X))
    worker <- function(j) {
      fit <- MuSiC::music.basic(as.numeric(Y[, j]), X, S, Sigma, iter.max = 1000, nu = 1e-4, eps = .01)
      z <- as.numeric(fit$p.weight); stopifnot(all(is.finite(z)), min(z) >= 0, sum(z) > 0)
      z / sum(z)
    }
    export <- c("X", "Y", "S", "Sigma"); columns <- colnames(X)
    kind <- "static_reference_basic_fraction_sensitivity"
  }
  cl <- parallel::makeCluster(cores)
  parallel::clusterCall(cl, function(paths) { invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8")); .libPaths(paths); NULL }, .libPaths())
  parallel::clusterExport(cl, c(export, "worker"), envir = environment())
  ans <- tryCatch(parallel::parLapply(cl, seq_len(ncol(bulk)), worker), finally = parallel::stopCluster(cl))
  pred <- do.call(rbind, ans); dimnames(pred) <- list(colnames(bulk), columns)
}
pred <- as.matrix(pred)
if (!all(colnames(bulk) %in% rownames(pred))) pred <- t(pred)
pred <- pred[colnames(bulk), setdiff(colnames(pred), "InfiltrationScore"), drop = FALSE]
stopifnot(nrow(pred) == 45L, all(is.finite(pred)), !anyDuplicated(colnames(pred)))
if (!grepl("score", kind)) stopifnot(min(pred) >= -1e-10, max(abs(rowSums(pred) - 1)) < 1e-5)
for (target.id in targets) {
  map <- mapping[mapping$method == task & mapping$target == target.id]
  stopifnot(all(parts(map$predicted_states) %in% colnames(pred)))
}
saved <- list(prediction = pred, method = task, output_kind = kind, details = details,
  protocol_md5 = protocol.hash, runner_md5 = md5(script), cores = cores,
  seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
saveRDS(saved, final)
write_mat(pred, paste0("predictions/", task, ".tsv"))
write_json(saved[setdiff(names(saved), "prediction")], paste0("predictions/", task, "_metadata.json"))
capture.output(sessionInfo(), file = file.path(out, "logs", paste0(task, "_sessionInfo.txt")))
cat("COMPLETE", task, round(saved$seconds, 2), "seconds\n")
