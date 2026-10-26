# Figure 3: equal-slice mean categorical concordance and rank-only RCTD contributions.
args<-commandArgs(TRUE);archive<-args[1];work<-args[2]
suppressPackageStartupMessages(library(ggplot2));suppressPackageStartupMessages(library(ggalign))
font<-systemfonts::match_fonts('Arial');stopifnot(grepl('arial',font$path,ignore.case=TRUE))
theme_set(theme_minimal(base_family='Arial',base_size=10)+theme(plot.margin=margin(2,3,2,2),plot.title=element_text(hjust=0),plot.subtitle=element_text(hjust=0),plot.title.position='plot',legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0)))
aliases<-c(Cellpose3='Cellpose cyto3',CellposeSAM='Cellpose cpsam',CellposeSAMv2='Cellpose cpsam v2',CellposeDINOsmall='Cellpose DINO small',CellposeDINO='Cellpose DINO',InstanSeg_prior_fluor='InstanSeg fluorescence',InstanSeg_HE='InstanSeg H&E nuclei',microSAM_LM_AIS='microSAM AIS',microSAM_LM_AMG='microSAM AMG',PathoSAM_HE_AIS='PathoSAM AIS',skimage_HE='scikit-image H&E',OpenCV_HE='OpenCV H&E',Zhang_initial_s2p5_nuclei='Watershed A: nuclei',Zhang_initial_s2p5_expand5='Watershed A: +5 µm',Zhang_initial_s3p5_nuclei='Watershed B: nuclei',Zhang_initial_s3p5_expand5='Watershed B: +5 µm',Zhang_revised_s2p5_nuclei='Watershed C: nuclei',Zhang_revised_s2p5_expand5='Watershed C: +5 µm',Zhang_revised_s3p5_nuclei='Watershed D: nuclei',Zhang_revised_s3p5_expand5='Watershed D: +5 µm')
prov<-file.path(archive,'provenance','figures','figure3')
a<-as.matrix(read.csv(file.path(prov,'CD1.agreement.csv'),row.names=1,check.names=FALSE))
b<-as.matrix(read.csv(file.path(prov,'CD2.agreement.csv'),row.names=1,check.names=FALSE))
stopifnot(identical(dimnames(a),dimnames(b)))
avg<-(a+b)/2;usable<-rowSums(is.finite(avg))>0;avg<-avg[usable,usable,drop=FALSE]
stopifnot(all(is.finite(avg)),max(abs(avg-t(avg)))<1e-12,all(diag(avg)==1))
write.csv(avg,file.path(work,'mean_agreement.csv'))
hc<-hclust(as.dist(1-avg),method='average');labs<-unname(aliases[rownames(avg)]);rownames(avg)<-colnames(avg)<-hc$labels<-labs
clean_dendro<-theme(axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank())
heat<-ggheatmap(avg)+coord_fixed(ratio=1)+
 scale_fill_gradientn(colours=c('#F7FBFF','#C6DBEF','#6BAED6','#2171B5','#08306B'),limits=c(0,1),labels=scales::label_percent(),name='Mean agreement',guide=guide_colourbar(direction='horizontal',barwidth=grid::unit(45,'mm'),barheight=grid::unit(3,'mm'),title.position='top'))+
 labs(x=NULL,y=NULL)+
 theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1,size=9),axis.text.y=element_text(size=9),plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=10),legend.position='bottom')+
 anno_top(size=grid::unit(6,'mm'))+align_dendro(distance=NULL,method=hc)+clean_dendro+
 labs(title='A  Primary-call concordance',subtitle='Mean across slices')+
 theme(plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=10),plot.title.position='plot')+
 anno_right(size=grid::unit(6,'mm'))+align_dendro(distance=NULL,method=rev(as.dendrogram(hc)))+clean_dendro+quad_active()
# The title is attached to the native top dendrogram annotation, above its tree.
z<-read.csv(file.path(work,'rank_by_slice.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(z)==360,length(unique(z$method))==20,setequal(unique(z$sample),c('CD1','CD2')))
methods<-unique(z$method);means<-do.call(rbind,lapply(methods,function(m){
 x<-z[z$method==m,];u<-x[x$sample=='CD1',];v<-x[x$sample=='CD2',];stopifnot(identical(u$rank,v$rank))
 data.frame(method=m,rank=u$rank,contribution=(u$contribution+v$contribution)/2,n=sum(c(u$n[1],v$n[1])))
}))
valid<-means[means$n>0,];stopifnot(all(valid$contribution>=0),all(abs(tapply(valid$contribution,valid$method,sum)-1)<1e-12))
first<-means[means$rank==1,];ordered<-first$method[order(-first$contribution,first$method,na.last=TRUE)]
method_ids<-setNames(sprintf('%02d',seq_along(aliases)),names(aliases))
bar_labels<-aliases
means$method_label<-factor(bar_labels[means$method],levels=bar_labels[ordered]);means$call<-factor(paste('Call',means$rank),levels=paste('Call',1:9))
write.csv(means,file.path(work,'mean_rank_contributions.csv'),row.names=FALSE)
missing<-means[means$rank==1 & means$n==0,];valid<-means[means$n>0,]
palette<-setNames(c('#009E73','#E69F00','#0072B2','#CC79A7','#D55E00','#56B4E9','#882255','#44AA99','#999999'),paste('Call',1:9))
bars<-ggplot(valid,aes(method_label,contribution,fill=call))+
 geom_col(position=position_stack(reverse=TRUE),width=.82,colour=NA)+
 geom_text(data=missing,aes(x=method_label,y=.5,label='No regions'),inherit.aes=FALSE,angle=90,size=3,family='Arial',colour='#666666')+
 scale_fill_manual(values=palette,drop=FALSE,name='Contribution rank')+
 scale_x_discrete(drop=FALSE)+scale_y_continuous(breaks=seq(0,1,.2),labels=scales::label_percent(),expand=expansion(mult=c(0,0)))+
 coord_cartesian(ylim=c(0,1))+
 labs(title='B  Ranked RCTD contributions',subtitle='Mean across slices; independent of cell identity',x=NULL,y='Mean relative contribution')+
 guides(fill=guide_legend(nrow=1,byrow=TRUE,title.position='top'))+
 theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=45,hjust=1,size=9),plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=10),legend.position='bottom',legend.text=element_text(size=7.5),legend.key.size=grid::unit(2.5,'mm'),legend.key.spacing.x=grid::unit(1,'mm'))
prior<-read.csv(file.path(prov,'compact_ranks','mean_rank_contributions.csv'),stringsAsFactors=FALSE)
stopifnot(identical(means$method,prior$method),identical(means$rank,prior$rank),isTRUE(all.equal(means$contribution,prior$contribution,tolerance=1e-14)))
prior_a<-as.matrix(read.csv(file.path(prov,'compact_ranks','mean_agreement.csv'),row.names=1,check.names=FALSE))
stopifnot(isTRUE(all.equal(unname(avg),unname(prior_a),tolerance=1e-14)))
# C pairs original RNA recovery (Figure 2 denominator) with the per-slice Call 1 mean.
rna<-read.csv(file.path(archive,'tables','rctd_full_major_summary.csv'),stringsAsFactors=FALSE)
fig2<-read.csv(file.path(archive,'tables/segmentation_metrics.csv'),stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM')
primary<-z[z$rank==1,c('method','sample','contribution','n')]
tradeoff<-merge(primary,rna[,c('method','sample','source_umi','assigned_umi','unassigned_umi','native_retained_regions')],by=c('method','sample'),all=TRUE,sort=FALSE)
stopifnot(nrow(tradeoff)==40,!anyDuplicated(paste(tradeoff$method,tradeoff$sample)),all(tradeoff$n==tradeoff$native_retained_regions),all(tradeoff$assigned_umi+tradeoff$unassigned_umi==tradeoff$source_umi))
idx<-match(paste(tradeoff$method,tradeoff$sample),paste(fig2$method,fig2$sample));stopifnot(!anyNA(idx))
tradeoff$assigned_fraction<-tradeoff$assigned_umi/tradeoff$source_umi
stopifnot(all(tradeoff$assigned_umi==fig2$assigned_umi[idx]),all(tradeoff$source_umi==fig2$source_umi[idx]),max(abs(100*tradeoff$assigned_fraction-fig2$assigned_umi_pct[idx]))<1e-9)
tradeoff$method_id<-unname(method_ids[tradeoff$method]);tradeoff$method_label<-unname(aliases[tradeoff$method])
tradeoff$included_in_c<-tradeoff$n>0 & is.finite(tradeoff$contribution)
by_slice<-tradeoff[tradeoff$included_in_c,];stopifnot(nrow(by_slice)==38,all(table(by_slice$sample)==19))
# Average the two slice-level metrics, never pool UMI denominators or polygons.
points<-do.call(rbind,lapply(split(by_slice,by_slice$method),function(d){
 stopifnot(nrow(d)==2,setequal(d$sample,c('CD1','CD2')))
 data.frame(method=d$method[1],method_id=d$method_id[1],method_label=d$method_label[1],assigned_fraction=mean(d$assigned_fraction),contribution=mean(d$contribution))
}))
stopifnot(nrow(points)==19,!anyDuplicated(points$method))
check_b<-first$contribution[match(points$method,first$method)];stopifnot(max(abs(points$contribution-check_b))<1e-14)
correlations<-data.frame(aggregation='arithmetic mean across slices',n=nrow(points),spearman_rho=cor(points$assigned_fraction,points$contribution,method='spearman'))
write.csv(tradeoff,file.path(work,'tradeoff_points.csv'),row.names=FALSE)
write.csv(points,file.path(work,'tradeoff_mean_points.csv'),row.names=FALSE)
write.csv(correlations,file.path(work,'tradeoff_correlations.csv'),row.names=FALSE)
c_methods<-names(aliases)[names(aliases)%in%points$method]
# Fixed qualitative colours selected for max-min CIEDE2000 separation and
# contrast >= 3 against white (minimum pair distance 16.10 vs 5.82 previously).
c_palette<-setNames(c('#222222','#37A61A','#D764D2','#D27B34','#009DC3','#5922AB','#8A1C07','#00522A','#786A39','#925A6A','#2769C5','#539E8A','#8C8DB7','#DF6F68','#777777','#11767C','#969529','#831464','#34466E'),c_methods)
stopifnot(length(c_palette)==19,!any(c_palette%in%palette))
scatter<-ggplot(points,aes(assigned_fraction,contribution,colour=method))+
 geom_point(size=3)+scale_colour_manual(values=c_palette,breaks=c_methods,labels=aliases[c_methods],name='Method')+
 scale_x_continuous(labels=scales::label_percent(),expand=expansion(mult=.12))+
 scale_y_continuous(labels=scales::label_percent(),expand=expansion(mult=.15))+
 labs(title='C  RNA recovery versus Call 1 dominance',subtitle=sprintf('Mean across slices | Spearman ρ = %.2f; n = %d methods',correlations$spearman_rho,correlations$n),x='Mean RNA assigned / original RNA',y='Mean Call 1 contribution')+
 guides(colour=guide_legend(ncol=1,override.aes=list(size=2.5)))+
 theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=10),legend.position='right',legend.text=element_text(size=7),legend.title=element_text(size=9),legend.key.height=grid::unit(2.5,'mm'),legend.key.width=grid::unit(3,'mm'),legend.key.spacing.y=grid::unit(.4,'mm'))
figure<-align_plots(heat,free_align(bars,axes='r'),free_align(scatter,axes='br'),design='AB\nAC',widths=c(1.1,1),heights=c(1.5,1),theme=theme(plot.caption=element_text(size=8,hjust=0),plot.margin=margin(2,2,2,2)))+
 layout_title(caption='All panels: arithmetic means of slice-level metrics. C: one point per method. RNA recovery includes all segmented regions;\nCall 1 dominance includes only RCTD-fitted polygons. Spearman correlation is descriptive, not an accuracy validation. microSAM AMG: no fitted regions.')
path<-file.path(work,'Figure3_RCTD_concordance_dominance.png')
ragg::agg_png(path,width=13,height=7.74,units='in',res=300,background='white');print(figure);invisible(dev.off())
jsonlite::write_json(list(status='RENDERED_PENDING_VISUAL_QA',path=path,md5=unname(tools::md5sum(path)),bytes=file.info(path)$size,dpi=300,font='Arial',font_path=font$path,heatmap='arithmetic mean of two per-slice categorical agreement matrices: (CD1 + CD2) / 2',bar='arithmetic mean of per-slice mean ranked relative coefficients',aggregation_unit='slice',patient_level_grouping=FALSE,numerical_values_unchanged=TRUE,call_identity_independent=TRUE,heatmap_square_constraint='native coord_fixed(ratio=1), height derived from width',row_dendrogram_side='right',column_dendrogram_side='top',no_thresholds=TRUE,no_clipping=TRUE,panel_c=list(points=38,methods_per_slice=19,x='assigned_umi/source_umi; exact Figure 2 metric',y='per-slice mean Call 1 relative contribution in fitted polygons',correlation='descriptive Spearman separately per slice; no pooled coefficient or inferential p-value',results=correlations)),file.path(work,'render.qa.json'),pretty=TRUE,auto_unbox=TRUE)
qa_path<-file.path(work,'render.qa.json');qa<-jsonlite::read_json(qa_path)
qa$panel_c<-list(points=19,one_point_per_method=TRUE,aggregation='arithmetic mean of the two slice-level fractions and Call 1 means; no patient grouping or pooling',x='mean assigned_umi/source_umi; exact Figure 2 metric per slice',y='mean of slice-level Call 1 relative contributions; identical to B',correlation='one descriptive Spearman across 19 method means; no inferential p-value',results=correlations)
qa$layout<-'native AB/AC; one heatmap, one rank bar chart, one scatter; data-adaptive scatter axes'
qa$appearance<-list(width_inches=13,height_inches=7.74,height_ratio_to_previous=.9,bar_axis_ids=FALSE,bar_legend_rows=1,heatmap_palette='original continuous Blues',heatmap_classes=NULL,heatmap_limits=c(0,1),scatter_palette='fixed max-min CIEDE2000 qualitative; contrast >=3 on white; distinct from B',scatter_palette_min_deltaE2000=16.09762,previous_palette_min_deltaE2000=5.816658,scatter_palette_values=c_palette,scatter_legend='right; method names',scatter_numeric_labels=FALSE,titles='all left aligned to native complete plot section',bar_right_space='independent of scatter legend; native free_align r')
jsonlite::write_json(qa,qa_path,pretty=TRUE,auto_unbox=TRUE)
cat('COMPACT_FIGURE3_RENDERED\n')
