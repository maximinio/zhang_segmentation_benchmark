# Validate native RCTD numeric output without clipping, replacing, or renormalizing it.
validate_native_rctd_weights<-function(weights,result){
 md<-as.data.frame(SummarizedExperiment::colData(result));tol<-sqrt(.Machine$double.eps)
 allowed_missing<-as.character(md$spot_class)=='reject'
 for(n in intersect(c('conv_all','conv_doublet'),names(md)))allowed_missing<-allowed_missing|(!is.na(md[[n]])&!md[[n]])
 allowed_missing[is.na(allowed_missing)]<-FALSE
 diagnostics<-lapply(names(weights),function(name){
  w<-weights[[name]]
  a<-as.matrix(w);stopifnot(identical(colnames(a),colnames(result)),!any(is.infinite(a)))
  badcols<-which(colSums(is.na(a))>0L);stopifnot(all(allowed_missing[badcols]))
  finite<-a[is.finite(a)];negative_columns<-which(colSums(a< -tol,na.rm=TRUE)>0L)
  # spacexr 1.4 create_spe_doublet stores unnormalised all_weights directly
  # as weights_full. Keep this auxiliary fit verbatim, NOT as probabilities.
  auxiliary<-identical(name,'weights_full')
  upper_columns<-which(colSums(a>1+tol,na.rm=TRUE)>0L)
  if(!auxiliary){stopifnot(all(allowed_missing[negative_columns]),all(allowed_missing[upper_columns]))}
  list(floating_point_tolerance=tol,negative_finite_entries=sum(finite<0),minimum_finite=if(length(finite))min(finite) else NULL,maximum_finite=if(length(finite))max(finite) else NULL,missing_entries=sum(is.na(a)),untrusted_native_columns=colnames(a)[sort(unique(c(badcols,negative_columns,upper_columns)))],native_auxiliary_raw_fit_not_probability_validated=auxiliary,primary_prediction_domain_checked=!auxiliary,native_values_preserved_without_clipping=TRUE)
 })
 names(diagnostics)<-names(weights)
 diagnostics
}
