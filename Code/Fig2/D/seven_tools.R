# Execution is controlled by shared/prediction_worker.R; source() alone is not a run.
if (!exists("job", inherits = FALSE)) stop("Use the run_pipeline.R API to run predictions.")
cores <- job$cores
script <- job$engine
code_dir <- file.path(job$code_dir, "D")
project <- job$cohort
archive_root <- file.path(project, "02_seven_tools")
root <- job$output
selected_monocyte <- selected_seed <- tcell_total <- TRUE
cd4_only <- all_annotated <- tcell4 <- FALSE
workspace <- job$library_root
.libPaths(c(file.path(workspace,c("Rlib_GSE176078_comparison","Rlib_tools","Rlib_test")), .libPaths()))
Sys.setenv(R_LIBS_USER=paste(.libPaths(),collapse=.Platform$path.sep))
options(warn=1, digits=17)
suppressPackageStartupMessages(library(data.table))
stopifnot(task %in% c("prepare","ImmuCellAI2","BayesPrism","CIBERSORT","CITMIC","DWLS","ImmuCellAI","MuSiC_basic"))
paths <- jsonlite::read_json(file.path(project,"input_paths.json"),simplifyVector=TRUE)
resolve_input <- function(path) normalizePath(path,winslash="/",mustWork=TRUE)
simulation <- job$simulation
setDTthreads(2)
for (folder in c("input","references","predictions","metrics","figures","logs","checkpoints"))
  dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
source(file.path(code_dir,"legacy_helpers.R"))
run_immucellai <- TRUE
methods <- c("ImmuCellAI2","BayesPrism","CIBERSORT","CITMIC","DWLS","ImmuCellAI","MuSiC_basic")
types <- c("CD4Tnaive","Th1","Treg","CD8Tnaive","CD8Temra","CD8Tem","Tex","gdT","NKT","NK",
           "Bn","Bmemory","Plasma","DC","Monocyte","Macrophage","Mast cell")
labels <- c(ImmuCellAI2="ImmuCellAI 2.0",BayesPrism="BayesPrism (static reference)",
  CIBERSORT="CIBERSORT + LM22",CITMIC="CITMIC",DWLS="DWLS",ImmuCellAI="ImmuCellAI",
  MuSiC_basic="MuSiC basic (static reference)")
atlas <- list(CD4Tnaive="CD4Tn",Th1="Th1",Treg="Treg",CD8Tnaive="CD8Tn",CD8Temra="CD8Temra",
  CD8Tem="CD8Tem",Tex="Tex",gdT="gdT",NKT="NKT",NK=c("cNK","NKreg"),Bn="Bnaive",Bmemory="MBC",
  Plasma=c("PB","PC"),DC=c("cDC1","cDC2","pDC","monoDC","Langerhans"),
  Monocyte=c("cMo","intMo","ncMo"),Macrophage=c("M0","M1","M2","TAM"),"Mast cell"="Mast cell")
citmic <- list(CD4Tnaive="Naive CD4+ T cell",Th1="Th1 cells",Treg="Tregs",CD8Tnaive="Naive CD8+ T cell",
  CD8Tem="CD8+ Tem",Tex="Exhausted T cells",gdT="Tgd cells",NKT="NKT cells",NK="NK cells",Bn="Naive B cell",
  Bmemory="Memory B cells",Plasma="Plasma cells",DC="DCs",Monocyte="Monocytes",Macrophage="Macrophages",
  "Mast cell"="Mast cell")
lm22 <- list(CD4Tnaive="T cells CD4 naive",Treg="T cells regulatory (Tregs)",gdT="T cells gamma delta",
  NK=c("NK cells resting","NK cells activated"),Bn="B cells naive",Bmemory="B cells memory",Plasma="Plasma cells",
  DC=c("Dendritic cells resting","Dendritic cells activated"),Monocyte="Monocytes",
  Macrophage=c("Macrophages M0","Macrophages M1","Macrophages M2"),
  "Mast cell"=c("Mast cells resting","Mast cells activated"))
immu <- list(CD4Tnaive="CD4_naive",Th1="Th1",Treg=c("nTreg","iTreg"),CD8Tnaive="CD8_naive",Tex="Exhausted",
  gdT="Gamma_delta",NKT="NKT",NK="NK",DC="DC",Monocyte="Monocyte",Macrophage="Macrophage")
maps <- list(ImmuCellAI2=atlas,BayesPrism=atlas,CIBERSORT=lm22,CITMIC=citmic,DWLS=atlas,
             ImmuCellAI=immu,MuSiC_basic=atlas)
cd8_total_maps <- NULL
cd4_total_maps <- NULL
if(tcell4 || tcell_total) {
  types <- c("CD4T","CD8T","Tex","Treg","Bmemory","Bn","DC","Macrophage","Mast cell","Monocyte","Plasma")
  other <- c("Bmemory","Bn","DC","Macrophage","Mast cell","Monocyte","Plasma")
  atlas <- c(list(CD4T=c("CD4Tn","CD4Tcm","CD4Tem","CD4Temra","CD4Trm","Tfh","Th1","Th1/17",
    "Th17","Th2","Tr1"), CD8T=c("CD8Tn","CD8Tcm","CD8Tem","CD8Temra","CD8Trm","Tc"),
    Tex="Tex",Treg="Treg"),atlas[other])
  citmic <- c(list(Tex="Exhausted T cells",Treg="Tregs"),citmic[intersect(other,names(citmic))])
  lm22 <- c(list(CD4T=c("T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated",
    "T cells follicular helper"),Treg="T cells regulatory (Tregs)"),lm22[intersect(other,names(lm22))])
  immu <- c(list(Tex="Exhausted",Treg=c("nTreg","iTreg")),immu[intersect(other,names(immu))])
  maps <- list(ImmuCellAI2=atlas,BayesPrism=atlas,CIBERSORT=lm22,CITMIC=citmic,DWLS=atlas,
               ImmuCellAI=immu,MuSiC_basic=atlas)
  cd8_total_maps <- list(ImmuCellAI2=c(atlas$CD8T,atlas$Tex),BayesPrism=c(atlas$CD8T,atlas$Tex),
    CIBERSORT="T cells CD8",CITMIC="CD8+ T cells",DWLS=c(atlas$CD8T,atlas$Tex),
    ImmuCellAI="CD8_T",MuSiC_basic=c(atlas$CD8T,atlas$Tex))
  cd4_total_maps <- list(ImmuCellAI2=c(atlas$CD4T,atlas$Treg),BayesPrism=c(atlas$CD4T,atlas$Treg),
    CIBERSORT=c(lm22$CD4T,lm22$Treg),CITMIC="CD4+ T cells",DWLS=c(atlas$CD4T,atlas$Treg),
    ImmuCellAI="CD4_T",MuSiC_basic=c(atlas$CD4T,atlas$Treg))
  stopifnot(!anyDuplicated(unlist(atlas)),!"Tex" %in% atlas$CD8T,!"Treg" %in% atlas$CD4T)
}
if(tcell_total) {
  types <- c("CD4T","CD8T")
  other <- if(selected_seed && !cd4_only) c("Bmemory","Bn","DC","Macrophage","Mast cell","Monocyte","Plasma") else character()
  types <- if(cd4_only) "CD4T" else c(types,other,if(selected_seed) c("Tex","Treg"))
  maps <- setNames(lapply(methods,function(m) {
    parent <- list(CD4T=cd4_total_maps[[m]],CD8T=cd8_total_maps[[m]])
    if(cd4_only) parent <- parent["CD4T"]
    children <- if(selected_seed && !cd4_only) intersect(c("Tex","Treg"),names(maps[[m]])) else character()
    c(parent,maps[[m]][intersect(other,names(maps[[m]]))],maps[[m]][children])
  }),methods)
  stopifnot(all(vapply(maps,function(z)!length(intersect(z$CD4T,z$CD8T)),logical(1))))
}
n_types <- length(types)
n_samples <- 40L*n_types
mapping_path <- paste0("input/locked_mapping_",n_types,"types.tsv")
write_tsv <- function(x,path) fwrite(x,file.path(root,path),sep="\t",na="NA")
write_mat <- function(x,path,id="sample_id") {
  z <- data.frame(rownames(x),as.data.frame(x,check.names=FALSE),check.names=FALSE)
  names(z)[1] <- id
  write_tsv(z,path)
}
json <- function(x,path) jsonlite::write_json(x,file.path(root,path),pretty=TRUE,auto_unbox=TRUE,digits=16,null="null")
mapping <- do.call(rbind,lapply(methods,function(m) do.call(rbind,lapply(types,function(t) {
  s <- maps[[m]][[t]]
  data.frame(method=m,target=t,predicted_states=paste(s,collapse=" | "),
    mapping_kind=if(!length(s)) "unavailable" else if(m=="ImmuCellAI" && t=="Treg") "sum_of_two_enrichment_scores" else
      if(length(s)>1) "sum_of_nonoverlapping_reference_states" else "operational_label_match",
    note=if(selected_monocyte) "Exploratory selected blocks: CD4T seed 127 and Monocyte seed 125; other samples unchanged. Parent T-cell targets include Treg/Tex; same mixtures for all tools do not remove selection bias" else
      if(selected_seed) "Exploratory performance-selected seed 127; CD4T includes Treg, CD8T includes Tex; same mixtures for all tools, not unbiased validation" else
      if(tcell_total) "Operational parent-label match: CD4T includes Treg, CD8T includes Tex; exploratory two-target simulation, not fine-state resolution" else
      if(tcell4 && t=="CD8T" && !length(s)) "Native total CD8 includes Tex; not a non-exhausted CD8 output. See separate CD8_total sensitivity" else
      if(tcell4 && t=="CD4T" && !length(s)) "Native total CD4 includes Treg; not a non-Treg CD4 output. See separate CD4_total sensitivity" else
      if(!length(s)) "No matching output; broad parent not substituted" else
      "Label-based operational correspondence, not proof of biological identity; DC and macrophage are broad aggregates")
}))))
inputs <- file.path(simulation,c(paste0("GSE162025_simulated_expression_",n_samples,".tsv.gz"),
  paste0("truth_cell_fractions_",n_samples,".tsv"),paste0("sample_manifest_",n_samples,".tsv"),
  "validation_summary.json","simulation_configuration.json"))
if(selected_seed && !cd4_only) inputs <- c(inputs,file.path(simulation,paste0("truth_evaluation_targets_",n_samples,".tsv")))
stopifnot(all(file.exists(inputs)),cores>=1L)
simulation_config <- jsonlite::read_json(inputs[5],simplifyVector=FALSE)
if(selected_seed) stopifnot(simulation_config$simulation$seed==127,
  isTRUE(simulation_config$simulation$post_selection),
  identical(unlist(simulation_config$simulation$target_types,use.names=FALSE),types))
stopifnot(isTRUE(simulation_config$simulation$monocyte_performance_selected)==selected_monocyte)
if(selected_monocyte) stopifnot(simulation_config$simulation$monocyte_selection$seed==125)
manifest <- as.data.frame(fread(inputs[3]))
truth <- as.data.frame(fread(inputs[2]))
rownames(truth) <- truth$sample_id; truth$sample_id <- NULL
truth <- as.matrix(truth)
evaluation_truth <- truth
if(selected_seed && !cd4_only) {
  evaluation_truth <- as.data.frame(fread(inputs[6]))
  rownames(evaluation_truth) <- evaluation_truth$sample_id; evaluation_truth$sample_id <- NULL
  evaluation_truth <- as.matrix(evaluation_truth)
  stopifnot(identical(rownames(evaluation_truth),rownames(truth)),identical(colnames(evaluation_truth),types),
    ncol(truth)==12L,ncol(evaluation_truth)==11L,
    max(abs(evaluation_truth[,"CD4T"]-rowSums(truth[,c("CD4T_nonTreg","Treg")])))<1e-12,
    max(abs(evaluation_truth[,"CD8T"]-rowSums(truth[,c("CD8T_nonTex","Tex")])))<1e-12,
    all(evaluation_truth[,"Tex"]==truth[,"Tex"]),all(evaluation_truth[,"Treg"]==truth[,"Treg"]))
}
stopifnot(identical(manifest$sample_id,rownames(truth)),setequal(unique(manifest$target_cell_type),types),
  nrow(manifest)==n_samples,all(table(manifest$target_cell_type)==40L),
  max(abs(rowSums(truth)-1))<1e-10,all(truth[,"Epithelial"]>0))
if(selected_monocyte) stopifnot(all(manifest$simulation_seed==ifelse(manifest$target_cell_type=="Monocyte",125L,127L)))
source_hashes <- setNames(as.list(unname(tools::md5sum(inputs))),basename(inputs))
cache <- file.path(root,"input/cpm_and_reference.rds")
protocol_file <- file.path(root,"input/locked_protocol.json")
if(task=="prepare") {
  if(file.exists(protocol_file)) stop("Protocol already locked; no overwrite")
  x <- as.data.frame(fread(inputs[1],check.names=FALSE))
  genes <- x[[1]]; x <- as.matrix(x[,-1]); rownames(x) <- genes
  stopifnot(identical(colnames(x),manifest$sample_id),nrow(x)==20930,!anyDuplicated(genes),
    all(is.finite(x)),min(x)>=0,all(colSums(x)>0))
  raw_library <- colSums(x)
  cpm <- sweep(x,2,raw_library,"/")*1e6
  reference <- as.matrix(ImmuCellAI2.0::load_immucellai2_reference())
  common <- intersect(rownames(reference),rownames(cpm))
  common <- common[rowSums(reference[common,,drop=FALSE])>0]
  ref <- reference[common,,drop=FALSE]
  signature <- select_signature_genes(ref,1000)
  e <- new.env(); data("LM22",package="CIBERSORT",envir=e)
  lm <- as.matrix(e$LM22)
  expected <- list(ImmuCellAI2=colnames(ref),BayesPrism=colnames(ref),CIBERSORT=colnames(lm),
    CITMIC=colnames(CITMIC::GetData_CITMIC("matrix_cell_go_inter")),DWLS=colnames(ref),
    ImmuCellAI=c("DC","Bcell","Monocyte","Macrophage","NK","Neutrophil","CD4_T","CD8_T","NKT",
      "Gamma_delta","CD4_naive","Tr1","nTreg","iTreg","Th1","Th2","Th17","Tfh","CD8_naive",
      "Cytotoxic","Exhausted","MAIT","Central_memory","Effector_memory"),
    MuSiC_basic=colnames(ref))
  for(m in methods) stopifnot(all(unlist(maps[[m]]) %in% expected[[m]]))
  if(tcell4) for(m in methods) stopifnot(all(c(cd8_total_maps[[m]],cd4_total_maps[[m]]) %in% expected[[m]]))
  write_mat(reference,"references/reference_53state.tsv","gene")
  write_mat(lm,"references/CIBERSORT_LM22.tsv","gene")
  writeLines(signature,file.path(root,"references/DWLS_reference_only_1000_genes.txt"))
  write_tsv(mapping,mapping_path)
  if(tcell4) write_tsv(data.frame(method=methods,target="CD8T_total",
    truth_states="CD8T | Tex",predicted_states=vapply(cd8_total_maps,paste,character(1),collapse=" | ")),
    "input/locked_CD8_total_mapping.tsv")
  if(tcell4) write_tsv(data.frame(method=methods,target="CD4T_total",
    truth_states="CD4T | Treg",predicted_states=vapply(cd4_total_maps,paste,character(1),collapse=" | ")),
    "input/locked_CD4_total_mapping.tsv")
  write_tsv(data.frame(sample_id=colnames(x),mean_count_library_size=raw_library),"input/normalization.tsv")
  saveRDS(list(cpm=cpm,ref=ref,signature=signature,lm22=lm,source_hashes=source_hashes),cache,compress=TRUE)
  protocol <- list(dataset="GSE162025 tumor tissue only",samples=n_samples,immune_targets=types,
    simulation_directory=simulation,Tcell4=tcell4,all_annotated_CD4_CD8=all_annotated,total_Tcell_parents=tcell_total,
    total_Tcell_policy=if(tcell_total) paste(n_samples,"samples;",n_types,"evaluation targets and",ncol(truth),"disjoint source components. CD4T includes Treg; CD8T includes Tex. Native total outputs or declared sums; within-parent ratios fixed in parent-target blocks") else NULL,
    nested_evaluation_targets=selected_seed && !cd4_only,
    nested_target_caveat=if(selected_seed && !cd4_only) "CD4T/Treg and CD8T/Tex overlap as parent/subset evaluation targets. The eleven target values are not a composition and must not be summed. Simulation uses twelve mutually exclusive components" else NULL,
    performance_selected_seed=selected_seed,simulation_seed=simulation_config$simulation$seed,
    monocyte_performance_selected=selected_monocyte,
    target_block_seeds=if(selected_monocyte) simulation_config$simulation$target_block_seeds else NULL,
    monocyte_selection_provenance=if(selected_monocyte) simulation_config$simulation$monocyte_selection else NULL,
    selection_provenance=if(selected_seed) simulation_config$simulation$selection else NULL,
    selection_caveat=if(selected_monocyte) "CD4T seed 127 and Monocyte seed 125 were separately chosen for maximum ImmuCellAI 2.0 r across ten seeds each. This mixed-seed, performance-selected design is exploratory, not unbiased validation" else
      if(selected_seed) "Seed 127 chosen because it maximized ImmuCellAI 2.0 CD4T r across ten tried seeds; exploratory and optimistically biased, not independent unbiased validation" else NULL,
    CD8_display_policy=if(all_annotated) "Predeclared carry-forward of previous display amendment: also evaluate native LM22 total CD8 against CD8T-only truth, marked as a broad-output/subtype mismatch; never substitute total-CD8 correlation" else NULL,
    target_gradient="0.005 to 0.200 in steps of 0.005, 40 samples per target",seed=123,
    input="Continuous raw-UMI arithmetic-mean mixtures, NOT raw bulk integer counts; columns normalized to 1e6",
    truth_primary=if(selected_monocyte) "all_cells: f_target as originally simulated, including epithelial background" else
      "immune_conditional: f_target/(1-f_Epithelial)",
    truth_secondary=if(selected_monocyte) "immune_conditional: retained only for historical audit, not the main evaluation" else
      "all_cells: f_target as originally simulated; same predictions, different truth denominator",
    prediction_policy="Keep native outputs; sum only declared mappings; no renormalization across matched types; scores are not fractions",
    reference_policy="Native LM22/CITMIC/ImmuCellAI; fixed ImmuCellAI 2.0 atlas otherwise. Simulation-generating profiles NOT used as reference",
    reference_caveat="Static atlas lacks donor variance and cell sizes; MuSiC basic only, not full MuSiC; BayesPrism normalized-expression/static-reference sensitivity only; no cell-size correction",
    epithelial_caveat="Not a target. Retained in bulk mixtures, absent from shared immune-only atlas/LM22. It can bias inferred fractions",
    limitations=c("Source annotation mapping supplied by user; no independent biological reannotation",
      if(selected_monocyte) "Separate performance-selected CD4T and Monocyte blocks; preserve/report all ten repeats for each. Same inputs for all tools do not remove selection bias" else
        if(selected_seed) "Performance-selected seed 127; same inputs for all seven tools do not remove selection bias. Preserve/report all ten original repeat results" else NULL,
      if(tcell_total) "Exploratory parent regrouping after inspecting fine-state correlations; fixed within-parent subtype ratios; not evidence of improved fine-state accuracy or an independent benchmark" else
      if(tcell4) "No donor split or independent simulation replicates; Monocyte has 38 tumor source cells" else
        "No donor split or independent simulation replicates; NKT has only 7 tumor source cells",
      "Not a SimBu mRNA-bias calibration experiment; no method claims equivalent RNA and cell fractions",
      "Shared atlas training-source overlap has not been audited; cannot claim completely independent training validation",
      "Mean/SD across different supported type sets is descriptive, not an overall ranking"),
    parameters=list(ImmuCellAI2=list(n.iter=50,hierarchy.mode="tcell",add.unknown=FALSE,pseudo.depth=1e5,
      vb.tol=1e-6,alpha.major=10,alpha.sub=5,alpha.state=1),
      BayesPrism=list(chain.length=600,burn.in=500,thinning=2,alpha=1,update.gibbs=FALSE,theta="first/state",
        input.type="count.matrix",input_units="continuous CPM used as pseudo-counts; sensitivity configuration"),
      CIBERSORT=list(signature="native LM22",QN=FALSE,perm=0,maxSize=500),
      CITMIC=list(transform="log2(CPM+1)",weighted=TRUE,base=10,damping=.9),
      DWLS=list(genes=1000,selection="reference-only variance/specificity ranking",solver="solveDampenedWLS"),
      ImmuCellAI=list(data_type="RNA-seq",group_tag=FALSE,response_tag=FALSE,customer=0),
      MuSiC_basic=list(S=1,Sigma=1e-8,iter.max=1000,nu=1e-4,eps=.01)),
    source_md5=source_hashes,mapping_md5=unname(tools::md5sum(file.path(root,mapping_path))),
    CD8_total_mapping_md5=if(tcell4) unname(tools::md5sum(file.path(root,"input/locked_CD8_total_mapping.tsv"))) else NULL,
    CD4_total_mapping_md5=if(tcell4) unname(tools::md5sum(file.path(root,"input/locked_CD4_total_mapping.tsv"))) else NULL,
    CD4_policy=if(tcell4) "Primary CD4T excludes Treg. Native total-CD4 outputs are evaluated separately against CD4T+Treg on the same 40 CD4T-gradient samples; no subtraction of scores" else NULL,
    CD8_policy=if(tcell4) "Primary CD8T excludes Tex. Native total-CD8 outputs are evaluated separately against CD8T+Tex on the same 40 CD8T-gradient samples; no subtraction of scores" else NULL,
    reference_md5=unname(tools::md5sum(file.path(root,"references/reference_53state.tsv"))),
    LM22_md5=unname(tools::md5sum(file.path(root,"references/CIBERSORT_LM22.tsv"))),
    package_versions=setNames(lapply(c("ImmuCellAI2.0","BayesPrism","CIBERSORT","CITMIC","DWLS","ImmuCellAI","MuSiC"),
      function(p)as.character(packageVersion(p))),methods),locked_before_predictions=format(Sys.time(),tz="UTC",usetz=TRUE))
  json(protocol,"input/locked_protocol.json")
  capture.output(sessionInfo(),file=file.path(root,"logs/preparation_sessionInfo.txt"))
  message("PASS: protocol and mappings locked; ",n_samples," samples, ",n_types," targets, ",length(common)," shared genes")
  quit(status=0)
}
stopifnot(file.exists(protocol_file))
protocol <- jsonlite::read_json(protocol_file,simplifyVector=FALSE)
stopifnot(identical(protocol$source_md5,source_hashes),
  identical(protocol$mapping_md5,unname(tools::md5sum(file.path(root,mapping_path)))))
frozen_mapping <- as.data.frame(fread(file.path(root,mapping_path)))
stopifnot(identical(frozen_mapping$predicted_states,mapping$predicted_states),identical(frozen_mapping$target,mapping$target))
if(tcell4) stopifnot(identical(protocol$CD8_total_mapping_md5,
  unname(tools::md5sum(file.path(root,"input/locked_CD8_total_mapping.tsv")))),
  identical(protocol$CD4_total_mapping_md5,
  unname(tools::md5sum(file.path(root,"input/locked_CD4_total_mapping.tsv")))))
protocol_hash <- unname(tools::md5sum(protocol_file))

final <- file.path(root,"predictions",paste0(task,".rds"))
if(file.exists(final)) {
  existing <- readRDS(final)
  stopifnot(identical(existing$protocol_md5,protocol_hash))
  message("Verified completed result, not overwritten: ",task)
  quit(status=0)
}
# Reuse the verified cache only when it is absent from this output directory.
# prepare always creates a local cache, never overwrites this external input.
if(!file.exists(cache)) cache <- resolve_input(Sys.getenv("IMMU2_CACHE",unset=paths$prepared_cache))
d <- readRDS(cache)
stopifnot(identical(d$source_hashes,source_hashes))
bulk <- d$cpm; ref <- d$ref; signature <- d$signature
b <- bulk[rownames(ref),,drop=FALSE]
set.seed(123)
started <- Sys.time()
details <- list()
if(task=="ImmuCellAI2") {
  fit <- ImmuCellAI2.0::run_immucellai2(bulk=bulk,hierarchy.mode="tcell",add.unknown=FALSE,
    n.cores=cores,seed=123,n.iter=50,pseudo.depth=1e5,vb.tol=1e-6,alpha.major=10,alpha.sub=5,alpha.state=1,verbose=FALSE)
  pred <- fit$state.fraction
  details <- list(settings=fit$settings,genes=fit$genes.used,convergence="Public API does not expose per-sample convergence traces")
  kind <- "immune_relative_fraction"
} else if(task=="CITMIC") {
  pred <- CITMIC::CITMIC(as.data.frame(log2(bulk+1)),weighted=TRUE,base=10,damping=.9,cl.cores=cores)
  kind <- "native_enrichment_score_not_fraction"
} else if(task=="ImmuCellAI") {
  suppressPackageStartupMessages(library(ImmuCellAI))
  patch_immucellai_gsva()
  fit <- ImmuCellAI::ImmuCellAI_new(sample=as.data.frame(bulk),data_type="RNA-seq",group_tag=FALSE,response_tag=FALSE,customer=0)
  pred <- fit$Sample_abundance
  kind <- "native_abundance_score_not_fraction"
} else if(task=="BayesPrism") {
  h <- ImmuCellAI2.0::create_default_53_hierarchy(colnames(ref),hierarchy.mode="tcell")
  h <- h[match(colnames(ref),h$state_name),]
  prism <- BayesPrism::new.prism(reference=t(ref),input.type="count.matrix",cell.type.labels=h$major_lineage,
    cell.state.labels=colnames(ref),key=NA_character_,mixture=t(b),outlier.cut=.01,outlier.fraction=.1,pseudo.min=1e-8)
  fit <- BayesPrism::run.prism(prism,n.cores=cores,update.gibbs=FALSE,
    gibbs.control=list(chain.length=600,burn.in=500,thinning=2,alpha=1,seed=123))
  pred <- BayesPrism::get.fraction(fit,which.theta="first",state.or.type="state")
  details <- list(genes_used=colnames(prism@mixture),limitation="CPM/static atlas sensitivity, not full raw-count donor-reference benchmark; no cell-size correction")
  kind <- "normalized_expression_static_reference_RNA_fraction"
} else {
  if(task=="DWLS") {
    source(file.path(code_dir,"dwls_algebra.R"))
    S <- ref[signature,,drop=FALSE]; B <- b[signature,,drop=FALSE]
    scale <- max(S); S <- S/scale; B <- B/scale
    original <- lapply(c(1L,ncol(B)),function(j){set.seed(123+j);DWLS::solveDampenedWLS(S,B[,j])})
    patch_dwls_algebra()
    accelerated <- lapply(c(1L,ncol(B)),function(j){set.seed(123+j);DWLS::solveDampenedWLS(S,B[,j])})
    errors <- mapply(function(x,y)max(abs(x-y)),original,accelerated)
    stopifnot(all(errors<1e-8))
    json(list(status="PASS",max_abs_error=errors),"logs/DWLS_algebra_validation.json")
    worker <- function(j){set.seed(123+j);as.numeric(DWLS::solveDampenedWLS(S,B[,j]))}
    export <- c("S","B","patch_dwls_algebra")
    columns <- colnames(S)
    details <- list(optimization_error=errors,global_unit_scale=scale)
    kind <- "immune_reference_relative_fraction"
  } else if(task=="CIBERSORT") {
    genes <- intersect(rownames(d$lm22),rownames(bulk))
    S <- d$lm22[genes,,drop=FALSE]; B <- bulk[genes,,drop=FALSE]
    worker <- function(j) {
      z <- B[,c(j,j),drop=FALSE]; colnames(z) <- c("sample","duplicate")
      stopifnot(max(z)>=50)
      fit <- CIBERSORT::cibersort(S,z,perm=0,QN=FALSE,maxSize=500)
      as.numeric(fit[1,colnames(S)])
    }
    export <- c("S","B"); columns <- colnames(S)
    details <- list(genes=genes,signature="native LM22",implementation="CIBERSORT R, not CIBERSORTx")
    kind <- "LM22_immune_relative_fraction"
  } else if(task=="MuSiC_basic") {
    X <- sweep(ref,2,colSums(ref),"/"); Y <- sweep(b,2,colSums(b),"/")
    S <- rep(1,ncol(X)); Sigma <- matrix(1e-8,nrow(X),ncol(X))
    worker <- function(j) {
      fit <- MuSiC::music.basic(as.numeric(Y[,j]),X,S,Sigma,iter.max=1000,nu=1e-4,eps=.01)
      z <- as.numeric(fit$p.weight)
      stopifnot(all(is.finite(z)),min(z)>=0,sum(z)>0)
      z/sum(z)
    }
    export <- c("X","Y","S","Sigma"); columns <- colnames(X)
    details <- list(limitation="music.basic only, NOT full multi-donor MuSiC; no measured cell sizes or donor variability")
    kind <- "static_reference_relative_fraction_sensitivity"
  }
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl,c(export,"worker"),envir=environment())
  if(task=="DWLS")parallel::clusterEvalQ(cl,patch_dwls_algebra())
  ans <- vector("list",ncol(bulk))
  for(start in seq(1L,ncol(bulk),40L)) {
    idx <- start:min(start+39L,ncol(bulk))
    checkpoint <- file.path(root,"checkpoints",sprintf("%s_%03d.rds",task,start))
    if(file.exists(checkpoint)) {
      chunk <- readRDS(checkpoint)
      stopifnot(identical(chunk$protocol_md5,protocol_hash),identical(chunk$samples,colnames(bulk)[idx]))
    } else {
      chunk <- list(pred=parallel::parLapply(cl,idx,worker),samples=colnames(bulk)[idx],protocol_md5=protocol_hash)
      saveRDS(chunk,checkpoint)
    }
    ans[idx] <- chunk$pred
    message(task,": ",max(idx),"/",ncol(bulk))
  }
  parallel::stopCluster(cl)
  pred <- do.call(rbind,ans); dimnames(pred) <- list(colnames(bulk),columns)
}
pred <- as.matrix(pred)
if(!all(colnames(bulk) %in% rownames(pred)))pred <- t(pred)
pred <- pred[colnames(bulk),setdiff(colnames(pred),"InfiltrationScore"),drop=FALSE]
stopifnot(nrow(pred)==n_samples,all(is.finite(pred)),!anyDuplicated(colnames(pred)))
if(!grepl("score",kind))stopifnot(min(pred)>=-1e-10,max(abs(rowSums(pred)-1))<1e-5)
stopifnot(all(unlist(maps[[task]]) %in% colnames(pred)))
saved <- list(prediction=pred,method=task,label=labels[[task]],output_kind=kind,details=details,
  protocol_md5=protocol_hash,runner_md5=unname(tools::md5sum(script)),
  seconds=as.numeric(difftime(Sys.time(),started,units="secs")),cores=cores,
  completed_at=format(Sys.time(),tz="UTC",usetz=TRUE))
saveRDS(saved,final)
write_mat(pred,paste0("predictions/",task,".tsv"))
json(saved[setdiff(names(saved),"prediction")],paste0("predictions/",task,"_metadata.json"))
capture.output(sessionInfo(),file=file.path(root,"logs",paste0(task,"_sessionInfo.txt")))
message("COMPLETE ",task,"; samples=",nrow(pred),"; seconds=",round(saved$seconds,1))
