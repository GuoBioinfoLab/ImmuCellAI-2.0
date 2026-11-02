# All CIBERSORT fits use the unmodified LM22 object bundled with the installed package.
read_benchmark_bulk <- function(path) {
  lines <- strsplit(readLines(path,n=2,warn=FALSE),'\t',fixed=TRUE)
  header <- lines[[1]]
  if(length(header)==length(lines[[2]])-1L) header<-c('gene',header)
  stopifnot(length(header)==length(lines[[2]]))
  z<-data.table::fread(path,skip=1,header=FALSE,col.names=header,check.names=FALSE)
  genes<-z[[1]]; m<-as.matrix(z[,-1,with=FALSE]); storage.mode(m)<-'double'; rownames(m)<-genes
  stopifnot(!anyDuplicated(genes),!anyDuplicated(colnames(m)),all(is.finite(m)),min(m)>=0)
  m
}
write_matrix <- function(m,path) {
  data.table::fwrite(data.frame(id=rownames(m),m,check.names=FALSE),path,sep='\t',na='NA')
}
lm22_mapper <- function(dataset) {
  if(dataset=='GSE176078') return(list(Bmemory='B cells memory',Bnaive='B cells naive',
    CD4Tnaive='T cells CD4 naive',M1='Macrophages M1',Monocyte='Monocytes',
    NK=c('NK cells resting','NK cells activated'),Tfh='T cells follicular helper',Treg='T cells regulatory (Tregs)'))
  if(dataset=='GSE146771') return(list(CD4Tmemory=c('T cells CD4 memory resting','T cells CD4 memory activated'),
    CD4Tn='T cells CD4 naive',macrophage=c('Macrophages M0','Macrophages M1','Macrophages M2'),
    mast_cell=c('Mast cells resting','Mast cells activated'),plasma='Plasma cells',Treg='T cells regulatory (Tregs)'))
  stopifnot(dataset=='GSE164522')
  m<-list(Bcells=c('B cells naive','B cells memory'),CD4Tfh='T cells follicular helper',
    CD4Tnaive='T cells CD4 naive',Monocyte='Monocytes',NK=c('NK cells resting','NK cells activated'),
    Treg='T cells regulatory (Tregs)',mast_cell=c('Mast cells resting','Mast cells activated'),plasma='Plasma cells')
  m[[paste0(intToUtf8(c(947,948)),'T')]]<-'T cells gamma delta'
  m
}
evaluate_lm22 <- function(pred,truth,dataset,out) {
  targets<-unique(truth$target_cell_type); mapping<-lm22_mapper(dataset)
  rows<-pairs<-crosswalk<-list()
  for(target in targets) {
    states<-mapping[[target]]
    if(length(states)) stopifnot(all(states %in% colnames(pred)))
    idx<-which(truth$target_cell_type==target)
    stopifnot(identical(truth$step[idx],1:40))
    x<-if(length(states)) rowSums(pred[idx,states,drop=FALSE]) else rep(NA_real_,40)
    y<-truth$target_fraction[idx]
    crosswalk[[target]]<-data.frame(dataset=dataset,target=target,supported=length(states)>0,
      predicted_states=paste(states,collapse=' | '),
      note=if(!length(states)) 'No matching LM22 resolution; do not substitute a broader parent or a differently defined state.' else
        if(target=='Bcells') 'Naive plus memory B cells; plasma cells remain separate.' else
        'Operational correspondence based on supplied labels; sibling outputs are summed without renormalizing.')
    pairs[[target]]<-data.frame(dataset=dataset,method='CIBERSORT_LM22',target=target,
      sample_id=truth$sample_id[idx],step=truth$step[idx],truth=y,prediction=x)
    for(stratum in c('Full','Low','Higher')) {
      keep<-switch(stratum,Full=rep(TRUE,40),Low=truth$step[idx]<=10,Higher=truth$step[idx]>10)
      xx<-x[keep]; yy<-y[keep]
      status<-if(!length(states)) 'unavailable' else if(max(xx)==min(xx)) 'constant_prediction' else 'valid'
      rows[[length(rows)+1L]]<-data.frame(dataset=dataset,method='CIBERSORT_LM22',target=target,stratum=stratum,
        n_designed=sum(keep),status=status,Pearson=if(status=='valid') cor(xx,yy) else NA_real_,
        predicted_states=paste(states,collapse=' | '),supported=length(states)>0)
    }
  }
  for(item in c('rows','pairs','crosswalk')) {
    name<-switch(item,rows='per_celltype_abundance_metrics.tsv',pairs='sample_level_truth_prediction.tsv',crosswalk='target_mapping.tsv')
    data.table::fwrite(do.call(rbind,get(item)),file.path(out,name),sep='\t',na='NA')
  }
}
run_lm22 <- function(base,dataset,cores=6L,output_dir) {
  stopifnot(dataset %in% c('GSE176078','GSE146771','GSE164522'),cores>=1L)
  root<-file.path(base,paste0('Fig2_',dataset)); out<-output_dir
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  bulk_file<-file.path(root,'input/bulk_original.txt')
  input_audit<-jsonlite::read_json(file.path(root,'audit/input_validation.json'))
  stopifnot(unname(tools::md5sum(bulk_file))==input_audit$source_md5)
  bulk<-read_benchmark_bulk(bulk_file)
  truth<-read.delim(file.path(root,'input/target_truth_manifest.tsv'),check.names=FALSE)
  stopifnot(identical(colnames(bulk),truth$sample_id),all(abs(truth$target_fraction-truth$step*.005)<1e-12))
  e<-new.env(); data('LM22',package='CIBERSORT',envir=e); ref<-as.matrix(e$LM22)
  stopifnot(nrow(ref)==547L,ncol(ref)==22L,!anyDuplicated(rownames(ref)),!anyDuplicated(colnames(ref)))
  write_matrix(ref,file.path(out,'LM22_reference.tsv'))
  genes<-intersect(rownames(ref),rownames(bulk)); S<-ref[genes,,drop=FALSE]; B<-bulk[genes,,drop=FALSE]
  stopifnot(length(genes)>100L,all(apply(B,2,max)>=50))
  writeLines(genes,file.path(out,'matched_LM22_genes.txt'))
  meta<-list(dataset=dataset,method='CIBERSORT R + LM22',package_version=as.character(packageVersion('CIBERSORT')),
    reference_source='data("LM22", package="CIBERSORT")',LM22_genes=547L,LM22_cell_types=22L,
    LM22_md5=unname(tools::md5sum(file.path(out,'LM22_reference.tsv'))),
    input_md5=unname(tools::md5sum(bulk_file)),matched_genes=length(genes),samples=ncol(bulk),
    QN=FALSE,perm=0,maxSize=500,cores=cores,seed=123,
    implementation_repository='https://github.com/Moonerss/CIBERSORT',official_CIBERSORTx=FALSE,
    signature_selection='LM22/input intersection only; no shared 53-state reference or top-1000 gene selection',
    sample_truth='Target fraction from 40-step manifest; other background fractions are not reconstructed.')
  desc<-packageDescription('CIBERSORT'); meta$installed_commit<-desc$RemoteSha
  prefix<-file.path(out,'CIBERSORT_LM22')
  if(file.exists(paste0(prefix,'_predictions.rds'))) {
    old<-jsonlite::read_json(paste0(prefix,'_metadata.json'))
    stopifnot(identical(old$input_md5,meta$input_md5),identical(old$LM22_md5,meta$LM22_md5),
      identical(old$package_version,meta$package_version),identical(old$QN,FALSE),old$perm==0)
    pred<-readRDS(paste0(prefix,'_predictions.rds'))
    stopifnot(identical(rownames(pred),truth$sample_id),identical(colnames(pred),colnames(ref)))
    evaluate_lm22(pred,truth,dataset,out)
    message(dataset,': revalidated existing corrected outputs'); return(invisible(meta))
  }
  worker<-function(j) {
    set.seed(123+j)
    z<-B[,c(j,j),drop=FALSE]; colnames(z)<-c('sample','duplicate')
    fit<-CIBERSORT::cibersort(S,z,perm=0,QN=FALSE,maxSize=500)
    stopifnot(max(abs(fit[1,colnames(S)]-fit[2,colnames(S)]))<1e-12)
    as.numeric(fit[1,colnames(S)])
  }
  environment(worker)<-.GlobalEnv
  cl<-parallel::makeCluster(cores); on.exit(parallel::stopCluster(cl),add=TRUE)
  parallel::clusterExport(cl,c('S','B'),envir=environment())
  started<-Sys.time(); ans<-vector('list',ncol(B))
  for(start in seq.int(1L,ncol(B),by=80L)) {
    idx<-start:min(ncol(B),start+79L)
    ans[idx]<-parallel::parLapply(cl,idx,worker)
    message(dataset,': CIBERSORT + LM22 ',max(idx),'/',ncol(B))
  }
  pred<-do.call(rbind,ans); dimnames(pred)<-list(colnames(B),colnames(ref))
  stopifnot(all(is.finite(pred)),min(pred)>=0,max(abs(rowSums(pred)-1))<1e-8)
  meta$seconds<-as.numeric(difftime(Sys.time(),started,units='secs'))
  meta$completed_at<-format(Sys.time(),tz='UTC',usetz=TRUE)
  if(dataset=='GSE176078') {
    prior<-readRDS(file.path(root,'predictions/CIBERSORT_full_predictions.rds'))
    stopifnot(identical(dimnames(pred),dimnames(prior)))
    meta$previous_LM22_max_abs_difference<-max(abs(pred-prior))
    stopifnot(meta$previous_LM22_max_abs_difference<1e-10)
  }
  saveRDS(pred,paste0(prefix,'_predictions.rds')); write_matrix(pred,paste0(prefix,'_predictions.tsv'))
  jsonlite::write_json(meta,paste0(prefix,'_metadata.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
  evaluate_lm22(pred,truth,dataset,out)
  capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
  message('COMPLETE ',dataset,' CIBERSORT + LM22; ',nrow(pred),' samples')
  invisible(meta)
}
