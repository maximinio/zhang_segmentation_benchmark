cfg<-jsonlite::fromJSON(commandArgs(TRUE)[1]);d<-cfg$stage
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(Matrix))
stopifnot(unname(tools::md5sum(cfg$source))==cfg$source_md5)
x<-readRDS(cfg$source);stopifnot(validObject(x),nrow(x)==18085,ncol(x)==988942,identical(Images(x),cfg$samples));m<-LayerData(x,assay='Spatial',layer='counts');stopifnot(all(m@x>=0),all(m@x==floor(m@x)))
inputs<-file.path(d,'inputs');dir.create(inputs,showWarnings=FALSE);qa<-list()
for(s in cfg$samples){
 im<-x[[s]];co<-GetTissueCoordinates(im,scale=NULL);co<-co[match(Cells(im),co$cell),];stopifnot(identical(co$cell,Cells(im)))
 b<-do.call(rbind,strsplit(sub('.*s_008um_','',co$cell),'_'));co$row<-as.integer(b[,1]);co$col<-as.integer(sub('-.*','',b[,2]));co$px<-co$x*im@scale.factors$hires;co$py<-co$y*im@scale.factors$hires
 bx<-lm(px~row+col,co);by<-lm(py~row+col,co);basis<-rbind(coef(bx)[c('row','col')],coef(by)[c('row','col')]);stopifnot(max(abs(c(resid(bx),resid(by))))<.25)
 co$inside_image<-co$px>=0&co$px<ncol(im@image)&co$py>=0&co$py<nrow(im@image);co$nUMI<-as.numeric(colSums(m[,co$cell,drop=FALSE]));co$nGenes<-as.integer(colSums(m[,co$cell,drop=FALSE]>0))
 p<-file.path(inputs,paste0(s,'.png'));stopifnot(!file.exists(p));png::writePNG(im@image,p);back<-png::readPNG(p);stopifnot(identical(dim(back),dim(im@image)),max(abs(back-im@image))<=1/255)
 write.csv(co,file.path(inputs,paste0(s,'.bins.csv')),row.names=FALSE)
 g<-list(sample=s,hires=im@scale.factors$hires,um_per_pixel=cfg$bin_um/mean(sqrt(colSums(basis^2))),basis=basis,shape=dim(im@image),lattice_fit_max_error=max(abs(c(resid(bx),resid(by)))),total_umi=sum(co$nUMI),inside_image_umi=sum(co$nUMI[co$inside_image]),total_bins=nrow(co),inside_image_bins=sum(co$inside_image),image_md5=unname(tools::md5sum(p)),native_image_class=class(im),coordinate_transform='original hires scale only; no coordinate fitting applied')
 stopifnot(round(3.5/g$um_per_pixel)==round(4.5/g$um_per_pixel));g$seed_setting_alias<-list(nominal_um=c(3.5,4.5),same_effective_min_distance_pixels=round(3.5/g$um_per_pixel))
 jsonlite::write_json(g,file.path(inputs,paste0(s,'.geometry.json')),auto_unbox=TRUE,pretty=TRUE);qa[[s]]<-g;cat('BENCHMARK_ORIGINAL_EXPORTED',s,nrow(co),sum(co$nUMI),'\n');flush.console()
}
jsonlite::write_json(list(valid=TRUE,source=cfg$source,source_md5=cfg$source_md5,source_catalog_id=271,genes=nrow(x),bins=ncol(x),sections=qa),file.path(inputs,'export.qa.json'),auto_unbox=TRUE,pretty=TRUE)
