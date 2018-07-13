#! /usr/bin/env Rscript

load('analysis/Rdata/deseq.sel.Rdata')

source('analysis/09-colors.R')

library(DESeq2)
library(pvclust)
library(dendextend)

makepdf <- TRUE

################################################################################
# Hierarchical clustering samples
################################################################################
set.seed(1234)

tmat <- as.matrix(assay(tform.sel))
samples.pv <- pvclust::pvclust(tmat,
                               method.dist = 'correlation',
                               use.cor = "pairwise.complete.obs",
                               method.hclust = 'average', 
                               nboot=1000
                               )

save(samples.pv, file='analysis/Rdata/pvclust_samples.Rdata')

if(makepdf) pdf('analysis/figures/sample_clustering.pdf', 
                paper='USr', width=11, height=7.5)

par(mar=c(7,3,7,1))

dorder <- samples.pv %>% as.dendrogram() %>% labels
dlabels <- colData(tform.sel)[dorder, 'display']
dcols <- celltype_colors[as.character(colData(tform.sel)[dorder,'celltype'])]

samples.pv %>% 
    as.dendrogram() %>% 
    dendextend::set('labels', dlabels) %>%
    dendextend::set('labels_cex', 0.7) %>%
    dendextend::set('labels_col', dcols) %>%
    hang.dendrogram() %>%    
    plot(main="Cluster dendrogram with AU/BP values")

samples.pv %>% text(cex=0.7) 

colored_dots(data.frame(
    "paired"=layout_colors[as.character(colData(tform.sel)[dorder,'layout'])],
    "center"=center_colors[as.character(colData(tform.sel)[dorder,'sequencing_center'])],
    "cancer"=karyotype_colors[as.character(colData(tform.sel)[dorder,'karyotype'])],
    celltype=dcols
), y_shift=-0.14)

if(makepdf) dev.off()

################################################################################
# Prune for celltypes dendrogram
################################################################################
dend.short <- samples.pv %>% 
    as.dendrogram() %>% 
    prune(dorder[colData(tform.sel)[dorder,'sequencing_center'] != 'CSHL'])

ds.labels <- as.character(colData(tform.sel)[labels(dend.short), 'celltype'])

dend.short <- dend.short %>%
    dendextend::set('labels', ds.labels) %>%
    dendextend::set('labels_col', celltype_colors[ds.labels])

if(makepdf) pdf('analysis/figures/celltype_clustering.pdf', 
                paper='USr', width=11, height=7.5)

    dend.short %>% plot

if(makepdf) dev.off()

celltypes_order <- ds.labels
save(dend.short, celltypes_order, file='analysis/Rdata/pvclust_celltypes.Rdata')
