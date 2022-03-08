#! /usr/bin/env Rscript

load('analysis/Rdata/deseq.sel.Rdata')
load('analysis/Rdata/pvclust_celltypes.Rdata')
source('analysis/09-colors.R')

library(tidyverse)
library(DESeq2)
library(cowplot)

makepdf <- TRUE

################################################################################
# PCA
################################################################################
ntop <- 500
rv <- matrixStats::rowVars(assay(tform.sel))
sel.gene <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(assay(tform.sel)[sel.gene, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

stopifnot(all(row.names(pca$x) == row.names(colData(tform.sel))))

pca.df <- colData(tform.sel) %>%
    data.frame %>%
    dplyr::mutate(celltype=factor(celltype, levels=celltypes_order)) %>%
    dplyr::mutate(sequencing=paste(sequencing_center, layout)) %>%
    dplyr::mutate(sequencing=factor(sequencing, levels=c('CSHL PAIRED', 'CALTECH SINGLE', 'CALTECH PAIRED'))) %>%
    cbind(pca$x)

p12 <- pca.df %>% 
    ggplot(aes_string(x = "PC1", y = "PC2", color = "celltype", shape='sequencing')) + 
    geom_point(size=3, alpha=0.75) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
    scale_color_manual(values=celltype_colors) +
    ggtitle('HERV expression PCA')

if(makepdf) pdf('analysis/figures/sample_pca_1v2.pdf', paper='USr', width=6, height=4)
p12
if(makepdf) dev.off()

p13 <- pca.df  %>%
    ggplot(aes_string(x = "PC1", y = "PC3", color = "celltype", shape='sequencing')) + 
    geom_point(size=3, alpha=0.75) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
    scale_color_manual(values=celltype_colors) +
    ggtitle('HERV expression PCA')

if(makepdf) pdf('analysis/figures/sample_pca_1v3.pdf', paper='USr', width=6, height=4)
p13
if(makepdf) dev.off()
