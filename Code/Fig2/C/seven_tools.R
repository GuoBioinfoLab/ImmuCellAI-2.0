# Execution is controlled by shared/prediction_worker.R; source() alone is not a run.
if (!exists("job", inherits = FALSE)) stop("Use the run_pipeline.R API to run predictions.")
method <- task
mode <- "full"; cores <- job$cores
methods <- c("ImmuCellAI2","BayesPrism","CIBERSORT","CITMIC","DWLS","ImmuCellAI","MuSiC_basic")
stopifnot(method %in% methods)
script <- job$engine
code_dir <- file.path(job$code_dir, "C")
cohort <- job$cohort
run_name <- "exploratory_Bmemory_cDC2_20261004"
selected_update <- TRUE
incremental <- FALSE
root <- job$output
isolated_lib <- file.path(job$library_root, "Rlib_GSE176078_comparison")
stopifnot(dir.exists(isolated_lib))
.libPaths(c(isolated_lib,.libPaths())); Sys.setenv(R_LIBS_USER=isolated_lib)
options(warn=1,digits=17)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ImmuCellAI2.0))
stopifnot(as.character(packageVersion("ImmuCellAI2.0"))=="0.1.8")
for(d in c("predictions","logs","references","checkpoints","metrics","plots"))
  dir.create(file.path(root,d),recursive=TRUE,showWarnings=FALSE)
design <- readRDS(job$design)
if(run_name %in% c("final_9targets_rebuilt", "exploratory_Bmemory_cDC2_20261004")) {
  stopifnot(length(design$targets)==9L,design$n_selected_samples==360L,
    "NK" %in% design$atoms,isTRUE(design$single_uniform_input),!isTRUE(design$uniform_background))
} else stopifnot(identical(design$seeds,123:142),length(design$targets)==12L,!"NK" %in% design$atoms)
profile_path <- if(!is.null(design$profile_file)) design$profile_file else file.path(root,"atomic_mean_count_profiles.rds")
stopifnot(unname(tools::md5sum(profile_path))==design$profile_md5,
  identical(unname(tools::md5sum(unlist(design$source_files))),unlist(design$source_md5,use.names=FALSE)))
profiles <- readRDS(profile_path)
reference_file <- file.path(cohort,"references/reference_53celltypesTPM20260518.txt")
stopifnot(unname(tools::md5sum(reference_file))=="5aa315cb9006fb95a3e5b0621bd1a923")
reference <- ImmuCellAI2.0::read_expression_matrix(reference_file)
stopifnot(isTRUE(all.equal(reference,ImmuCellAI2.0::load_immucellai2_reference())))
write_tsv <- function(x,path) data.table::fwrite(x,path,sep="\t",na="NA")
write_mat <- function(x,path,id="sample_id") {
  z <- data.frame(rownames(x),x,check.names=FALSE); names(z)[1] <- id; write_tsv(z,path)
}
parameters <- list(hierarchy.mode="tcell",add.unknown=FALSE,n.cores=2L,seed=123L,
  n.iter=50L,pseudo.depth=1e5,vb.tol=1e-6,alpha.major=10,alpha.sub=5,alpha.state=1,verbose=FALSE)
source(file.path(code_dir,"legacy_helpers.R"))
run_immucellai <- TRUE
set.seed(123)
data_dir <- job$selected_dir
input <- file.path(data_dir,"simulated_expression.tsv.gz")
manifest <- read.delim(file.path(data_dir,"sample_manifest.tsv"),check.names=FALSE)
selection <- jsonlite::read_json(file.path(data_dir,"selection.json"),simplifyVector=TRUE)
stopifnot(unname(tools::md5sum(input))==selection$expression_md5)
bulk <- as.matrix(read.delim(gzfile(input),row.names=1,check.names=FALSE))
stopifnot(ncol(bulk)==design$n_selected_samples,identical(colnames(bulk),manifest$sample_id),all(is.finite(bulk)),min(bulk)>=0)
bulk <- sweep(bulk,2,colSums(bulk),"/")*1e6
if(mode=="smoke") bulk <- bulk[,c(1,40,ncol(bulk)-39,ncol(bulk)),drop=FALSE]
common <- intersect(rownames(reference),rownames(bulk))
common <- common[rowSums(reference[common,,drop=FALSE])>0]
ref <- as.matrix(reference[common,,drop=FALSE]); b <- bulk[common,,drop=FALSE]
signature <- select_signature_genes(ref,1000)
stopifnot(ncol(ref)==53L,length(signature)==1000L)

pkg <- switch(method, ImmuCellAI2 = "ImmuCellAI2.0", MuSiC_basic = "MuSiC", method)
metadata <- list(method = method, mode = mode, seed = 123, cores = cores,
  package = pkg, version = as.character(packageVersion(pkg)), samples = colnames(bulk),
  input_scale = "Verified raw-count mean-profile mixtures, normalized to CPM before inference",
  selection = design$selection,
  simulation_design_md5 = unname(tools::md5sum(job$design)),
  library_path = isolated_lib,
  bulk_md5 = unname(tools::md5sum(input)),
  runner_md5 = unname(tools::md5sum(script)),
  reference_md5 = unname(tools::md5sum(reference_file)), common_genes = length(common),
  reference_policy = "Native signatures when supplied; otherwise frozen ImmuCellAI 2.0 53-state atlas",
  historical_reference_source_overlap = "not newly audited")
if (selected_update) {
  metadata$post_selection <- TRUE
  metadata$limitation_selection <- design$caveat
}
if (incremental) {
  stopifnot(mode == "full")
  previous_prefix <- file.path(design$previous_run, "predictions", paste0(method, "_full"))
  previous_pred <- readRDS(paste0(previous_prefix, "_predictions.rds"))
  previous_meta <- jsonlite::read_json(paste0(previous_prefix, "_metadata.json"), simplifyVector = TRUE)
  previous_input <- file.path(design$previous_run, "selected_360", "simulated_expression.tsv.gz")
  stopifnot(previous_meta$bulk_md5 == unname(tools::md5sum(previous_input)),
            previous_meta$version == metadata$version)
  previous_bulk <- as.matrix(read.delim(gzfile(previous_input), row.names = 1, check.names = FALSE))
  previous_bulk <- sweep(previous_bulk, 2, colSums(previous_bulk), "/") * 1e6
  changed_indices <- which(manifest$target_cell_type %in% design$replaced_targets)
  retained_indices <- setdiff(seq_len(ncol(bulk)), changed_indices)
  retained_ids <- colnames(bulk)[retained_indices]
  replay_index <- retained_indices[40L]
  replay_id <- colnames(bulk)[replay_index]
  work_indices <- sort(c(changed_indices, replay_index))
  stopifnot(length(changed_indices) == 80L, length(retained_indices) == 280L,
    identical(rownames(previous_bulk), rownames(bulk)),
    identical(bulk[, retained_ids], previous_bulk[, retained_ids]),
    identical(match(retained_ids, colnames(previous_bulk)), retained_indices))
  metadata$incremental <- list(previous_run = design$previous_run,
    previous_prediction_md5 = unname(tools::md5sum(paste0(previous_prefix, "_predictions.rds"))),
    previous_metadata_md5 = unname(tools::md5sum(paste0(previous_prefix, "_metadata.json"))),
    new_samples = colnames(bulk)[changed_indices], replay_sample = replay_id,
    unchanged_samples = retained_ids, retained_input_exactly_identical = TRUE,
    note = "Unchanged-sample cache reuse requires identical preprocessing and a native inference replay. New selected blocks are recomputed for each tool.")
}
write_mat(reference, file.path(root, "references/reference_53state.tsv"), "gene")
writeLines(signature, file.path(root, "references/DWLS_reference_selected_1000_genes.txt"))
prefix <- file.path(root, "predictions", paste(method, mode, sep = "_"))
if (file.exists(paste0(prefix, "_predictions.rds"))) {
  previous <- jsonlite::read_json(paste0(prefix,"_metadata.json"),simplifyVector=TRUE)
  stopifnot(previous$bulk_md5==metadata$bulk_md5,previous$runner_md5==metadata$runner_md5)
  message("Verified completed method: ",method); quit(status=0)
}
started <- Sys.time()
if (method == "ImmuCellAI2") {
  metadata$settings <- "Default 53-state reference and marker genes; hierarchy.mode=tcell; add.unknown=FALSE; seed=123"
  fit <- ImmuCellAI2.0::run_immucellai2(bulk = bulk, hierarchy.mode = "tcell",
    add.unknown = FALSE, n.cores = cores, seed = 123, n.iter=50,
    pseudo.depth=1e5, vb.tol=1e-6, verbose = FALSE)
  pred <- fit$state.fraction
  metadata$fit_settings <- fit$settings
  metadata$genes_used <- fit$genes.used
  metadata$output_kind <- "relative_fraction"
} else if (method == "DWLS") {
  source(file.path(code_dir, "dwls_algebra.R"))
  S <- ref[signature, , drop = FALSE]
  B <- b[signature, , drop = FALSE]
  scale <- max(S)
  S <- S / scale
  B <- B / scale
  # Validate the optimization against native DWLS on this very input before use.
  validation_file <- file.path(root, "logs", paste0("DWLS_", mode, "_algebra_validation.json"))
  original <- lapply(c(1L, ncol(B)), function(j) {
    set.seed(123 + j)
    DWLS::solveDampenedWLS(S, B[, j])
  })
  patch_dwls_algebra()
  accelerated <- lapply(c(1L, ncol(B)), function(j) {
    set.seed(123 + j)
    DWLS::solveDampenedWLS(S, B[, j])
  })
  errors <- mapply(function(x, y) max(abs(x - y)), original, accelerated)
  stopifnot(all(errors < 1e-8))
  jsonlite::write_json(list(status = "PASS", samples = colnames(B)[c(1L, ncol(B))],
    max_abs_error = errors), validation_file, pretty = TRUE, auto_unbox = TRUE)
  metadata$settings <- "Native solveDampenedWLS; 1000 reference-only selected genes; joint global unit scaling; verified algebraic optimization, original bootstrap/damping/QP retained"
  metadata$common_unit_scale <- scale
  metadata$output_kind <- "relative_fraction"
  metadata$optimization_validation_max_error <- errors
  worker <- function(j) {
    set.seed(123 + j)
    as.numeric(DWLS::solveDampenedWLS(S, B[, j]))
  }
  export <- c("S", "B", "patch_dwls_algebra")
} else if (method == "CIBERSORT") {
  e <- new.env()
  data("LM22", package = "CIBERSORT", envir = e)
  native_ref <- as.matrix(e$LM22)
  genes <- intersect(rownames(native_ref), rownames(bulk))
  S <- native_ref[genes, , drop = FALSE]
  B <- bulk[genes, , drop = FALSE]
  write_mat(native_ref, file.path(root, "references/CIBERSORT_LM22.tsv"), "gene")
  metadata$reference_md5 <- unname(tools::md5sum(file.path(root, "references/CIBERSORT_LM22.tsv")))
  metadata$common_genes <- length(genes)
  metadata$settings <- "CIBERSORT R package, NOT CIBERSORTx; native LM22; QN=FALSE for RNA-seq; perm=0; maxSize=500"
  metadata$output_kind <- "relative_fraction"
  # A two-column call preserves matrix dimensions in the native implementation.
  worker <- function(j) {
    z <- B[, c(j, j), drop = FALSE]
    colnames(z) <- c("sample", "duplicate")
    stopifnot(max(z) >= 50)
    fit <- CIBERSORT::cibersort(S, z, perm = 0, QN = FALSE, maxSize = 500)
    as.numeric(fit[1, colnames(S)])
  }
  export <- c("S", "B")
} else if (method == "MuSiC_basic") {
  X <- sweep(ref, 2, colSums(ref), "/")
  Y <- sweep(b, 2, colSums(b), "/")
  S <- rep(1, ncol(X))
  Sigma <- matrix(1e-8, nrow(X), ncol(X))
  metadata$settings <- "Exploratory MuSiC::music.basic only, NOT full multi-donor MuSiC; static atlas; S=1; Sigma=1e-8; iter.max=1000; nu=1e-4; eps=0.01; no surrogate fallback"
  metadata$limitation <- "No donor variance or cell-size information is available in the static atlas. Normalized-expression input. Not a full-method benchmark."
  metadata$output_kind <- "relative_fraction_static_reference_sensitivity"
  worker <- function(j) {
    fit <- MuSiC::music.basic(as.numeric(Y[, j]), X, S, Sigma,
      iter.max = 1000, nu = 1e-4, eps = 0.01)
    z <- as.numeric(fit$p.weight)
    stopifnot(all(is.finite(z)), min(z) >= 0, sum(z) > 0)
    z / sum(z)
  }
  export <- c("X", "Y", "S", "Sigma")
} else if (method == "CITMIC") {
  metadata$reference_md5 <- NA_character_
  metadata$settings <- "Native CITMIC signatures; log2(supplied linear expression+1); weighted=TRUE; base=10; damping=0.9"
  metadata$output_kind <- "native_enrichment_score_not_cell_fraction"
  pred <- CITMIC::CITMIC(as.data.frame(log2(bulk + 1)), weighted = TRUE,
    base = 10, damping = 0.90, cl.cores = cores)
} else if (method == "ImmuCellAI") {
  suppressPackageStartupMessages(library(ImmuCellAI))
  patch_immucellai_gsva()
  metadata$reference_md5 <- NA_character_
  metadata$settings <- "Native signatures; ImmuCellAI_new RNA-seq defaults; GSVA API compatibility adapter; original score compensation/rounding; no across-label renormalization"
  metadata$output_kind <- "native_abundance_score_not_cell_fraction"
  fit <- ImmuCellAI::ImmuCellAI_new(sample = as.data.frame(bulk), data_type = "RNA-seq",
    group_tag = FALSE, response_tag = FALSE, customer = 0)
  pred <- fit$Sample_abundance
} else if (method == "BayesPrism") {
  metadata$settings <- "BayesPrism native Gibbs inference; static 53-state atlas and normalized linear expression (count.matrix API); first state-level theta; chain.length=600; burn.in=500; thinning=2; alpha=1; seed=123; update.gibbs=FALSE"
  metadata$limitation <- "Normalized-expression/static-reference sensitivity configuration, not raw-count single-cell-reference analysis; no cell-size correction"
  metadata$output_kind <- "relative_fraction_normalized_expression_sensitivity"
  h <- ImmuCellAI2.0::create_default_53_hierarchy(colnames(ref), hierarchy.mode = "tcell")
  h <- h[match(colnames(ref), h$state_name), ]
  prism <- BayesPrism::new.prism(reference = t(ref), input.type = "count.matrix",
    cell.type.labels = h$major_lineage, cell.state.labels = colnames(ref), key = NA_character_,
    mixture = t(b), outlier.cut = .01, outlier.fraction = .1, pseudo.min = 1e-8)
  if (incremental) {
    previous_prism <- BayesPrism::new.prism(reference = t(ref), input.type = "count.matrix",
      cell.type.labels = h$major_lineage, cell.state.labels = colnames(ref), key = NA_character_,
      mixture = t(previous_bulk[common, , drop = FALSE]),
      outlier.cut = .01, outlier.fraction = .1, pseudo.min = 1e-8)
    stopifnot(identical(prism@phi_cellState@phi, previous_prism@phi_cellState@phi),
      identical(prism@phi_cellType@phi, previous_prism@phi_cellType@phi),
      identical(prism@mixture[retained_ids, ], previous_prism@mixture[retained_ids, ]))
    # Native initial Gibbs uses the same seed separately for each sample; no shared reference update.
    native_body <- paste(deparse(body(BayesPrism:::run.gibbs.refPhi.ini)), collapse = "\n")
    stopifnot(grepl("set.seed(seed)", native_body, fixed = TRUE))
    metadata$incremental$full_batch_preprocessing_identical <- TRUE
    metadata$incremental$native_independent_sampler <- native_body
    jsonlite::write_json(list(passed = TRUE, phi_identical = TRUE, retained_mixtures_identical = TRUE,
      gene_count = ncol(prism@mixture), new_samples = 80L, replay_samples = 1L),
      file.path(root, "logs/BayesPrism_incremental_preflight.json"), auto_unbox = TRUE, pretty = TRUE)
    if (Sys.getenv("IMMU2_PREFLIGHT_ONLY") == "1") quit(status = 0)
    prism@mixture <- prism@mixture[colnames(bulk)[work_indices], , drop = FALSE]
  }
  fit <- BayesPrism::run.prism(prism, n.cores = cores, update.gibbs = FALSE,
    gibbs.control = list(chain.length = 600, burn.in = 500, thinning = 2, alpha = 1, seed = 123))
  pred <- BayesPrism::get.fraction(fit, which.theta = "first", state.or.type = "state")
  if (incremental) {
    error <- max(abs(pred[replay_id, ] - previous_pred[replay_id, colnames(pred)]))
    stopifnot(is.finite(error), error < 1e-10)
    merged <- matrix(NA_real_, ncol(bulk), ncol(pred), dimnames = list(colnames(bulk), colnames(pred)))
    merged[retained_ids, ] <- previous_pred[retained_ids, colnames(pred)]
    merged[rownames(pred), ] <- pred
    pred <- merged
    metadata$incremental$native_replay_max_error <- error
  }
} else stop("Unknown method: ", method)

if (method %in% c("DWLS", "CIBERSORT", "MuSiC_basic")) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, c(export, "worker"), envir = environment())
  if (method == "DWLS") parallel::clusterEvalQ(cl, patch_dwls_algebra())
  ans <- vector("list", ncol(bulk))
  if (incremental) {
    output_columns <- if (method == "CIBERSORT") colnames(S) else colnames(ref)
    stopifnot(all(output_columns %in% colnames(previous_pred)),
      previous_meta$reference_md5 == metadata$reference_md5)
    for (j in retained_indices) ans[[j]] <- as.numeric(previous_pred[colnames(bulk)[j], output_columns])
  } else work_indices <- seq_len(ncol(bulk))
  checkpoint_dir <- file.path(root, "checkpoints", paste(method, mode, sep = "_"))
  dir.create(checkpoint_dir, showWarnings = FALSE)
  for (start in seq(1L, length(work_indices), by = 40L)) {
    idx <- work_indices[start:min(start + 39L, length(work_indices))]
    path <- file.path(checkpoint_dir, sprintf("batch_%03d.rds", start))
    if (file.exists(path)) {
      cached <- readRDS(path)
      stopifnot(identical(cached$bulk_md5, metadata$bulk_md5),
        identical(cached$runner_md5, metadata$runner_md5),
        identical(cached$reference_md5, metadata$reference_md5),
        identical(cached$samples, colnames(bulk)[idx]))
      chunk <- cached$results
    } else {
      chunk <- parallel::parLapply(cl, idx, worker)
      saveRDS(list(samples = colnames(bulk)[idx], bulk_md5 = metadata$bulk_md5,
        reference_md5 = metadata$reference_md5, runner_md5 = metadata$runner_md5,
        results = chunk), path)
    }
    ans[idx] <- chunk
    message(method, ": checkpoint ", min(start + 39L, length(work_indices)), "/", length(work_indices))
  }
  parallel::stopCluster(cl)
  pred <- do.call(rbind, ans)
  dimnames(pred) <- list(colnames(bulk), if (method == "CIBERSORT") colnames(S) else colnames(ref))
  if (incremental) {
    error <- max(abs(pred[replay_id, ] - previous_pred[replay_id, colnames(pred)]))
    stopifnot(is.finite(error), error < 1e-8)
    metadata$incremental$native_replay_max_error <- error
    metadata$incremental$preprocessing <- "Reference-only signature selection and scaling; sample-wise normalization; unchanged inputs verified; DWLS retains original full-matrix sample-index seeds."
  }
}
pred <- as.matrix(pred)
if (!all(colnames(bulk) %in% rownames(pred)) && all(colnames(bulk) %in% colnames(pred))) pred <- t(pred)
pred <- pred[colnames(bulk), setdiff(colnames(pred), "InfiltrationScore"), drop = FALSE]
stopifnot(nrow(pred) == ncol(bulk), all(is.finite(pred)), !anyDuplicated(colnames(pred)))
if (grepl("^relative_fraction", metadata$output_kind))
  stopifnot(min(pred) >= -1e-10, max(abs(rowSums(pred) - 1)) < 1e-5)
metadata$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
metadata$completed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
metadata$columns <- colnames(pred)
metadata$prediction_range <- range(pred)
saveRDS(pred, paste0(prefix, "_predictions.rds"))
write_mat(pred, paste0(prefix, "_predictions.tsv"))
jsonlite::write_json(metadata, paste0(prefix, "_metadata.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", digits = 16)
capture.output(sessionInfo(), file = paste0(prefix, "_sessionInfo.txt"))
message("COMPLETE ", method, ": ", nrow(pred), " samples; seconds=", metadata$seconds)
