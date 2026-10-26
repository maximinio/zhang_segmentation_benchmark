# Three-stage workflow, using unchanged native Seurat polygons and native RCTD.
# Prepare exact accepted-mask geometry first with zhang_figure1_prepare.py,
# then run: Rscript zhang_figure1_render.R <benchmark archive> <work directory>
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
archive <- normalizePath(args[1], mustWork = TRUE)
work <- normalizePath(args[2], mustWork = TRUE)
out <- file.path(archive, 'figures', 'main')
dir.create(out, showWarnings = FALSE)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SpatialExperiment))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(sf))
sample <- 'CD2'
method <- 'Zhang_revised_s3p5_expand5'
job <- paste(method, sample, sep = '__')
roi <- jsonlite::fromJSON(file.path(archive, 'provenance/visual_fields.json'))[[sample]][1, ]
palette <- c(Myeloid='#D55E00', CD4T='#0072B2', CD8T='#56B4E9', ILC='#CC79A7',
             B_Plasma='#882255', Epithelial='#E69F00', Mesenchymal='#009E73',
             Endothelial='#44AA99', Neural='#AA4499', 'Not fitted'='#B3B3B3',
             'Tied maximum'='#222222')
font <- systemfonts::match_fonts('Arial')$path
stopifnot(grepl('arial', font, ignore.case=TRUE))
cache <- file.path(work, 'native_field.rds')
if (!file.exists(cache)) {
  cfg <- jsonlite::fromJSON(file.path(archive, 'provenance/zhang_benchmark_config.json'))
  stopifnot(unname(tools::md5sum(cfg$source)) == cfg$source_md5)
  geometry_qa <- jsonlite::fromJSON(file.path(work,'field_geometry.qa.json'))
  stopifnot(geometry_qa$method==method, geometry_qa$sample==sample,
            geometry_qa$exact_pixel_area, geometry_qa$all_parts_and_holes_preserved)
  cat('READ_ORIGINAL_SEURAT_AND_ACCEPTED_D_RESULT\n'); flush.console()
  original <- readRDS(cfg$source)
  geom <- st_read(file.path(work,'field.polygons.geojson'),quiet=TRUE)
  st_crs(geom) <- NA
  stopifnot(all(st_is_valid(geom)),
            all(as.numeric(st_area(geom))==geom$area_hires_px2))
  tab <- read.csv(file.path(archive, 'tables/rctd_full_major_summary.csv'), stringsAsFactors=FALSE)
  row <- tab[tab$job == job, ]
  stopifnot(nrow(row)==1L, unname(tools::md5sum(row$path))==row$md5)
  b <- readRDS(row$path)
  stopifnot(validObject(b$native_result), b$mode=='full', b$labels=='author major_cluster',
            b$reference_assay=='RNAraw')
  md <- as.data.frame(colData(b$input_regions))
  idx <- match(geom$label,md$source_label)
  stopifnot(!anyNA(idx),!anyDuplicated(md$source_label))
  ids <- md$region_id[idx]
  im <- original[[sample]]
  fullres <- geom
  st_geometry(fullres) <- st_geometry(geom)/im@scale.factors$hires
  spatial <- sp::spChFIDs(as(fullres,'Spatial'),ids)
  seg <- as(as(spatial,'SpatialPolygons'),'Segmentation')
  stopifnot(identical(seg@polygons,spatial@polygons),identical(comment(spatial),'TRUE'))
  comment(seg) <- comment(spatial)
  stopifnot(max(abs(as.numeric(st_area(st_as_sf(seg)))-as.numeric(st_area(fullres))))<1e-4)
  xy <- spatialCoords(b$input_regions)[ids,,drop=FALSE]
  im@boundaries <- list(segmentations=seg,
    centroids=CreateCentroids(data.frame(x=xy[,1],y=xy[,2],cell=ids)))
  im@assay <- 'Spatial.Polygons'; DefaultBoundary(im) <- 'segmentations'
  small <- CreateSeuratObject(assay(b$input_regions,'counts')[,ids,drop=FALSE],
    assay='Spatial.Polygons',meta.data=md[ids,,drop=FALSE],min.cells=0,min.features=0)
  small[[sample]] <- im
  stopifnot(identical(small@images[[sample]]@boundaries$segmentations@polygons,seg@polygons),
            identical(small[[sample]]@image,original[[sample]]@image))
  comment(small@images[[sample]]@boundaries$segmentations) <- comment(seg)
  w <- as.matrix(assay(b$native_result, 'weights'))
  groups <- names(palette)[1:9]
  stopifnot(setequal(rownames(w), groups), all(is.finite(w)))
  w <- w[groups, , drop=FALSE]
  maximum <- apply(w, 2, max)
  top <- max.col(t(w), ties.method='first')
  calls <- rownames(w)[top]
  calls[colSums(w==rep(maximum, each=nrow(w)))>1L] <- 'Tied maximum'
  md$display <- 'Not fitted'
  md$display[match(colnames(w), rownames(md))] <- calls
  idx <- match(small$source_label, md$source_label)
  stopifnot(!anyNA(idx), !anyDuplicated(md$source_label))
  small$workflow_major <- factor(md$display[idx], levels=names(palette))
  small$workflow_constant <- 'Region'
  geometry <- jsonlite::fromJSON(file.path(archive, paste0('provenance/',sample,'.geometry.json')))
  info <- list(source_md5=cfg$source_md5, native_bundle_md5=row$md5,
               source_mask_md5=geometry_qa$mask_md5,
               sample=sample, method=method, roi=roi, regions=length(ids),
               exact_native_polygons=TRUE, unchanged_native_coefficients=TRUE,
               um_per_pixel=geometry$um_per_pixel)
  saveRDS(list(object=small, info=info), cache, compress=FALSE)
  rm(original, b, geom, spatial, w, md); gc()
  cat('NATIVE_FIELD_CACHED',length(ids),'\n'); flush.console()
}
cached <- readRDS(cache); small <- cached$object; info <- cached$info
stopifnot(validObject(small), info$sample==sample, info$method==method,
          identical(as.numeric(unlist(info$roi[c('x','y','size')])),
                    as.numeric(unlist(roi[c('x','y','size')]))))
make_plot <- function(group, cols, alpha, stroke, outline=NULL) {
  p <- SpatialDimPlot(small, images=sample, group.by=group,
        image.scale='hires', plot_segmentations=TRUE, crop=TRUE,
        image.alpha=1, alpha=alpha, stroke=stroke, stroke.alpha=1,
        cols=cols, combine=FALSE)[[1]]
  polygon_layer <- which(vapply(p$layers, function(l) inherits(l$geom,'GeomPolygon'), logical(1)))
  stopifnot(length(polygon_layer)==1L)
  if (!is.null(outline)) p$layers[[polygon_layer]]$aes_params$colour <- outline
  p$layers[[polygon_layer]]$show.legend <- identical(group, 'workflow_major')
  scale_length <- 100/info$um_per_pixel
  p + coord_fixed(xlim=c(roi$x,roi$x+roi$size),
                  ylim=c(roi$y,roi$y+roi$size), expand=FALSE) +
    annotate('segment', x=roi$x+15, xend=roi$x+15+scale_length,
             y=roi$y+roi$size-22, yend=roi$y+roi$size-22,
             colour='white', linewidth=2) +
    annotate('text', x=roi$x+15, y=roi$y+roi$size-38,
             label='100 µm', hjust=0, colour='white', size=3.2, family='Arial') +
    theme_void(base_family='Arial') +
    theme(plot.title=element_text(size=15, face='bold', hjust=0, margin=margin(b=5)),
          plot.subtitle=element_text(size=10, hjust=0, margin=margin(b=8)),
          plot.caption=element_text(size=10.5, hjust=0, lineheight=1.13, margin=margin(t=11)),
          plot.margin=margin(5,4,5,4), legend.position='bottom')
}
a <- make_plot('workflow_constant', c(Region='transparent'), 0, 0) +
  labs(title='A  Input', subtitle='Original H&E',
       caption='Deposited 8 µm bin-level Seurat\nCounts, coordinates and H&E') +
  guides(fill='none', colour='none', alpha='none')
b <- make_plot('workflow_constant', c(Region='transparent'), 0, .22, '#28CF73') +
  labs(title='B  Segmentation', subtitle='Accepted region boundaries',
       caption='20 segmentation configurations\nWhole-bin RNA assignment by bin centre') +
  guides(fill='none', colour='none', alpha='none')
c <- make_plot('workflow_major', palette, .70, 0) +
  labs(title='C  RCTD annotation', subtitle='Dominant major-group component', fill=NULL,
       caption='Full-mode RCTD; scIBD reference\n(Nie et al., 2023)\nLargest native coefficient') +
  scale_fill_manual(values=palette, drop=FALSE,
                    labels=function(x) gsub('_',' / ',x)) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE, override.aes=list(alpha=1, linewidth=0))) +
  theme(legend.text=element_text(size=9), legend.key.size=unit(3,'mm'),
        legend.margin=margin(5,0,0,0), legend.spacing.x=unit(2,'mm'))
arrow_plot <- ggplot() +
  annotate('segment',x=.08,xend=.92,y=.5,yend=.5,linewidth=.75,colour='#4B5563',
           arrow=grid::arrow(length=unit(2.5,'mm'),type='closed')) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE) +
  theme_void(base_family='Arial') + theme(plot.margin=margin(0,0,0,0))
p <- wrap_plots(a,arrow_plot,b,arrow_plot,c,nrow=1,widths=c(1,.12,1,.12,1),guides='collect') +
  plot_annotation(caption=paste0('Starting from deposited 8 µm bins, without Space Ranger reprocessing.\n',
                   'Same CD2 field in all panels; illustrative configuration: H&E watershed D (+5 µm).'),
                  theme=theme(text=element_text(family='Arial'),
                    plot.caption=element_text(size=10,hjust=0,lineheight=1.1,margin=margin(t=9)),
                    plot.margin=margin(8,10,8,10))) & theme(legend.position='bottom')
target <- file.path(out,'Figure1_segmentation_RCTD_workflow.png')
ggsave(target,p,width=13,height=5.8,dpi=300,device=ragg::agg_png,bg='white',limitsize=FALSE)
info$path <- target; info$md5 <- unname(tools::md5sum(target)); info$dpi <- 300
info$font <- 'Arial'; info$font_path <- font; info$native_plot_api <- 'Seurat::SpatialDimPlot'
info$manual_coordinate_transformation <- FALSE
info$no_new_segmentation_or_RCTD <- TRUE; info$visual_inspection <- 'pending'
jsonlite::write_json(info,file.path(archive,'provenance','figures','figure1','render.qa.json'),pretty=TRUE,auto_unbox=TRUE)
cat('FIGURE1_RENDERED',target,'\n')
