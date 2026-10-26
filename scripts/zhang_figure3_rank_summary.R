# Summarise ordered RCTD contributions, independent of cell-type identity.
args<-commandArgs(TRUE);archive<-args[1];work<-args[2]
suppressPackageStartupMessages(library(SpatialExperiment))
dir.create(work,recursive=TRUE,showWarnings=FALSE)
tab<-read.csv(file.path(archive,'tables','rctd_full_major_summary.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(tab)==40L,all(as.logical(tab$readback_exact)),all(as.logical(tab$counts_conserved_exact)))
out<-list();audit<-list()
for(i in seq_len(nrow(tab))){
 r<-tab[i,];cat('READ',r$job,'\n');flush.console()
 receipt<-jsonlite::fromJSON(file.path(archive,'provenance','rctd',paste0('fit_',r$job,'.exit.json')))
 stopifnot(receipt$exit_code==0,file.info(r$path)$size==r$bytes,unname(tools::md5sum(r$path))==r$md5)
 b<-readRDS(r$path);stopifnot(b$mode=='full',b$labels=='author major_cluster',b$reference_assay=='RNAraw')
 v<-rep(NA_real_,9);n<-0L;neg<-0L
 if(!is.null(b$native_result)){
  w<-as.matrix(assay(b$native_result,'weights'));n<-ncol(w)
  stopifnot(nrow(w)==9,n==r$native_retained_regions,all(is.finite(w)),all(colSums(w)>0))
  neg<-sum(w<0);relative<-sweep(w,2,colSums(w),'/')
  ranked<-apply(relative,2,sort,decreasing=TRUE);v<-rowMeans(ranked)
  stopifnot(all(diff(v)<=0),abs(sum(v)-1)<1e-12,all(v>=0))
  audit[[r$job]]<-list(n=n,source_md5=r$md5,negative_native_entries_preserved=neg,min_native=min(w),display_sum_error=max(abs(colSums(relative)-1)))
 }else stopifnot(r$input_regions==0,r$native_retained_regions==0)
 out[[i]]<-data.frame(method=r$method,sample=r$sample,rank=1:9,contribution=v,n=n)
 rm(b);if(exists('w'))rm(w,relative,ranked);gc()
 cat('RANK_MEANS',r$job,n,paste(signif(v,4),collapse=','),'\n');flush.console()
}
z<-do.call(rbind,out);write.csv(z,file.path(work,'rank_by_slice.csv'),row.names=FALSE)
jsonlite::write_json(list(status='ALL40_RANK_SUMMARIES_VERIFIED',audits=audit,total_fitted_polygons=sum(tab$native_retained_regions),identity_independent=TRUE,normalization='native coefficient / unmodified column sum, no clipping',aggregation='mean per rank within each slice; equal mean across CD1 and CD2 in renderer'),file.path(work,'rank.qa.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('ALL40_RANK_SUMMARIES_VERIFIED\n')
