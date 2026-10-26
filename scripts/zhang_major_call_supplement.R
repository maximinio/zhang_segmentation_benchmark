# Extract display-only argmax calls from unchanged native full-mode RCTD fits.
args<-commandArgs(TRUE);archive<-args[1];work<-args[2]
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SpatialExperiment))
tab<-read.csv(file.path(archive,'tables','rctd_full_major_summary.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(tab)==40L,all(tab$readback_exact),all(tab$counts_conserved_exact))
groups<-c('Myeloid','CD4T','CD8T','ILC','B_Plasma','Epithelial','Mesenchymal','Endothelial','Neural')
dir.create(file.path(work,'calls'),showWarnings=FALSE)
for(i in seq_len(nrow(tab))){
 row<-tab[i,];cat('READ_NATIVE',row$job,'\n');flush.console()
 receipt<-jsonlite::fromJSON(file.path(archive,'provenance','rctd',paste0('fit_',row$job,'.exit.json')))
 stopifnot(receipt$exit_code==0L,file.info(row$path)$size==row$bytes,unname(tools::md5sum(row$path))==row$md5)
 b<-readRDS(row$path);input<-b$input_regions
 stopifnot(validObject(input),b$mode=='full',b$labels=='author major_cluster',b$reference_assay=='RNAraw')
 md<-as.data.frame(colData(input));stopifnot(identical(rownames(md),colnames(input)),nrow(md)==row$input_regions,!anyDuplicated(md$source_label))
 md$display_major<-rep('Not fitted',nrow(md));md$max_native_coefficient<-rep(NA_real_,nrow(md))
 if(!is.null(b$native_result)){
  stopifnot(validObject(b$native_result),identical(assayNames(b$native_result),'weights'))
  w<-as.matrix(assay(b$native_result,'weights'));stopifnot(setequal(rownames(w),groups),all(is.finite(w)),ncol(w)==row$native_retained_regions)
  w<-w[groups,,drop=FALSE];idx<-match(colnames(w),rownames(md));stopifnot(!anyNA(idx),!anyDuplicated(idx))
  maximum<-apply(w,2,max);top<-max.col(t(w),ties.method='first');tie<-colSums(w==rep(maximum,each=nrow(w)))>1L
  stopifnot(identical(as.numeric(w[cbind(top,seq_len(ncol(w)))]),as.numeric(maximum)))
  call<-rownames(w)[top];call[tie]<-'Tied maximum'
  md$display_major[idx]<-call;md$max_native_coefficient[idx]<-maximum
  stopifnot(all(md$rctd_status[idx]=='native_full_fit'),sum(md$display_major!='Not fitted')==ncol(w))
 }else stopifnot(nrow(md)==0L,row$native_retained_regions==0L)
 write.csv(md[,c('region_id','source_label','display_major','max_native_coefficient','rctd_status')],file.path(work,'calls',paste0(row$job,'.csv')),row.names=FALSE,na='')
 cat('ARGMAX_DISPLAY_EXACT',row$job,nrow(md),sum(md$display_major!='Not fitted'),sum(md$display_major=='Tied maximum'),'\n');flush.console()
 rm(b,input,md);if(exists('w'))rm(w);gc()
}
source<-file.path(archive,'inputs','Zhang2025_CD_fibrostenosis_VisiumHD_GSE302569_counts_only.rds')
stopifnot(unname(tools::md5sum(source))=='762bde2b546922d698290bc6f9a93018')
original<-readRDS(source)
for(s in c('CD1','CD2')){
 p<-file.path(work,paste0(s,'.png'));png::writePNG(original[[s]]@image,p)
 geom<-jsonlite::fromJSON(file.path(archive,'provenance',paste0(s,'.geometry.json')))
 stopifnot(unname(tools::md5sum(p))==geom$image_md5)
 cat('ORIGINAL_HE_PNG_MD5_EXACT',s,'\n');flush.console()
}
cat('ALL_40_NATIVE_ARGMAX_EXPORTS_VERIFIED\n')
