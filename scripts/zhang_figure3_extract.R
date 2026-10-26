args<-commandArgs(TRUE);archive<-args[1];work<-args[2]
suppressPackageStartupMessages(library(SpatialExperiment))
tab<-read.csv(file.path(archive,'tables','rctd_full_major_summary.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(tab)==40L,all(tab$readback_exact),all(tab$counts_conserved_exact))
groups<-c('Myeloid','CD4T','CD8T','ILC','B_Plasma','Epithelial','Mesenchymal','Endothelial','Neural')
dir.create(file.path(work,'profiles'),showWarnings=FALSE)
audit<-list()
for(i in seq_len(nrow(tab))){
 row<-tab[i,];cat('READ',row$job,'\n');flush.console()
 receipt<-jsonlite::fromJSON(file.path(archive,'provenance','rctd',paste0('fit_',row$job,'.exit.json')))
 stopifnot(receipt$exit_code==0L,file.info(row$path)$size==row$bytes,unname(tools::md5sum(row$path))==row$md5)
 b<-readRDS(row$path);md<-as.data.frame(colData(b$input_regions))
 stopifnot(b$mode=='full',b$labels=='author major_cluster',b$reference_assay=='RNAraw',nrow(md)==row$input_regions)
 calls<-data.frame(source_label=md$source_label,call=integer(nrow(md)))
 out<-list(job=row$job,method=row$method,sample=row$sample,w=matrix(numeric(),9,0,dimnames=list(groups,NULL)))
 if(!is.null(b$native_result)){
  w<-as.matrix(assay(b$native_result,'weights'))[groups,,drop=FALSE]
  stopifnot(all(is.finite(w)),ncol(w)==row$native_retained_regions,all(colSums(w)>0))
  ix<-match(colnames(w),rownames(md));stopifnot(!anyNA(ix))
  top<-max.col(t(w),ties.method='first');ties<-colSums(w==rep(apply(w,2,max),each=9))>1
  stopifnot(!any(ties));calls$call[ix]<-top
  # Display-only division by the native column sum. No clipping, threshold or alteration of the source.
  p<-sweep(w,2,colSums(w),'/');dominance<-apply(p,2,max)
  ord<-order(-dominance,colnames(w),method='radix')
  out$w<-p[,ord,drop=FALSE];out$source_label<-md$source_label[ix][ord]
  audit[[row$job]]<-list(n=ncol(w),native_min=min(w),native_max=max(w),negative_entries=sum(w<0),relative_sum_error=max(abs(colSums(p)-1)),dominance_range=range(dominance),input_md5=row$md5)
  cat('PROFILE',row$job,ncol(w),'MIN',min(w),'NEG',sum(w<0),'\n');flush.console()
 }
 saveRDS(out,file.path(work,'profiles',paste0(row$job,'.rds')),compress=FALSE)
 write.csv(calls,file.path(work,'profiles',paste0(row$job,'.csv')),row.names=FALSE)
 rm(b,md,out,calls);if(exists('w'))rm(w,p);gc()
}
jsonlite::write_json(audit,file.path(work,'extract.qa.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('ALL40_EXTRACTED\n')
