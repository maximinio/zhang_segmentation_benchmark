args<-commandArgs(TRUE);cfg<-jsonlite::fromJSON(args[1]);job<-args[2];d<-cfg$stage
suppressPackageStartupMessages(library(spacexr));suppressPackageStartupMessages(library(SpatialExperiment));suppressPackageStartupMessages(library(Matrix))
stopifnot(cfg$rctd_mode=='full',cfg$reference_labels=='major_cluster',cfg$reference_assay=='RNAraw');set.seed(cfg$seed)
q<-jsonlite::fromJSON(file.path(d,'counts',paste0(job,'.qa.json')));stopifnot(q$counts_readback_exact,q$per_gene_counts_conserved_exact,unname(tools::md5sum(q$path))==q$md5)
input<-readRDS(q$path);stopifnot(validObject(input));ids<-colnames(input);result<-NULL;native_filter<-NULL
if(length(ids)>0){
 rq<-jsonlite::fromJSON(file.path(d,'reference.qa.json'));stopifnot(rq$raw_counts_exact,unname(tools::md5sum(file.path(d,'reference.rds')))==rq$md5)
 reference<-readRDS(file.path(d,'reference.rds'))
 cat('FULL_MAJOR_NATIVE_PREPROCESS',job,ncol(input),'\n');flush.console()
 data<-createRctd(input,reference,cell_type_col=cfg$reference_labels,require_int=TRUE)
 cat('FULL_MAJOR_NATIVE_FIT',job,'\n');flush.console()
 result<-runRctd(data,rctd_mode='full',max_cores=as.integer(cfg$cores))
 stopifnot(validObject(result),all(colnames(result)%in%ids),identical(assayNames(result),'weights'))
 w<-as.matrix(assay(result,'weights'));stopifnot(all(is.finite(w)))
 tol<-1e-6;sums<-colSums(w)
 numeric_qa<-list(native_min=min(w),native_max=max(w),negative_entries=sum(w< -tol),above_one_entries=sum(w>1+tol),column_sum_max_error=max(abs(sums-1)),domain_ok=all(w>= -tol)&all(w<=1+tol)&all(abs(sums-1)<=tol),native_values_preserved_without_clipping=TRUE)
}else{
 numeric_qa<-list(native_min=NULL,native_max=NULL,negative_entries=NULL,above_one_entries=NULL,column_sum_max_error=NULL,domain_ok=NA,native_values_preserved_without_clipping=TRUE)
}
status<-rep('excluded_by_native_RCTD_filters',length(ids));names(status)<-ids
if(!is.null(result))status[colnames(result)]<-'native_full_fit'
colData(input)$rctd_status<-status
bundle<-list(native_result=result,input_regions=input,reference_audit=jsonlite::fromJSON(file.path(d,'reference.qa.json')),mode='full',labels='author major_cluster',reference_assay='RNAraw',numeric_qa=numeric_qa,package_version=as.character(packageVersion('spacexr')),interpretation='Native full-mode weights and metadata preserved verbatim. Exclusions are native RCTD preprocessing, not new QC/cell calling. No singlet/doublet confidence or rejection labels are invented for full mode. Deconvolution is not validation of boundaries.')
out<-file.path(cfg$destination,'results',job);dir.create(out,recursive=TRUE,showWarnings=FALSE)
target<-file.path(out,'native_full_major.rds');pending<-file.path(out,'native_full_major.pending.rds');stopifnot(!file.exists(target),!file.exists(pending))
saveRDS(bundle,pending,compress=FALSE);back<-readRDS(pending);stopifnot(identical(back,bundle),validObject(back$input_regions));if(!is.null(result))stopifnot(validObject(back$native_result))
stopifnot(file.rename(pending,target))
meta<-metadata(input);summary<-list(job=job,method=as.character(colData(input)$method[1]),sample=as.character(colData(input)$sample_id[1]),status=if(length(ids))'NATIVE_FULL_FIT_READBACK_VERIFIED' else 'ZERO_REGIONS_NO_RCTD_FIT',input_regions=length(ids),zero_RNA_regions=q$zero_RNA_regions,source_umi=q$source_umi,assigned_umi=q$assigned_umi,unassigned_umi=q$unassigned_umi,native_retained_regions=if(is.null(result))0L else ncol(result),native_excluded_regions=length(ids)-if(is.null(result))0L else ncol(result),retained_assigned_umi=if(is.null(result))0 else sum(colData(input)[colnames(result),'nUMI']),cell_types=if(is.null(result))NULL else rownames(result),numeric_qa=numeric_qa,readback_exact=TRUE,counts_conserved_exact=TRUE,path=target,md5=unname(tools::md5sum(target)),bytes=file.info(target)$size)
jsonlite::write_json(summary,file.path(out,'qa.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
cat('FULL_MAJOR_RESULT_READBACK_PASS',job,summary$native_retained_regions,'\n');flush.console()
