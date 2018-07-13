#! /usr/bin/env Rscript
load('analysis/Rdata/loaded.Rdata')

library(tidyverse)
library(DESeq2)


do.parallel <- FALSE
if(do.parallel) {
    library(BiocParallel)
    register(MulticoreParam(8))
}

makepdf <- TRUE


################################################################################
# Best counts - DESeq2 (PolyA extract only)
################################################################################
# Criteria for samples
extracts <- c("longPolyA")
sfilt <- samples$rnaextract %in% extracts

# Criteria for genes
mincount <- 5
minsamp  <- 2
gfilt <- rowSums(counts.best[,sfilt] >= mincount) >= minsamp
gfilt <- gfilt & (annot.herv$chrom != 'chrY')

dds.sel <- DESeqDataSetFromMatrix(counts.best[gfilt, sfilt], samples[sfilt, ], ~celltype)
sizeFactors(dds.sel) <- calculateSizeFactors(metrics[sfilt,]$mapped_frags.bt2)
dds.sel <- DESeq(dds.sel, parallel=do.parallel)
tform.sel <- varianceStabilizingTransformation(dds.sel, blind=FALSE)
tform.sel.0 <- varianceStabilizingTransformation(dds.sel, blind=TRUE)

################################################################################
# Best counts - Hierarchical clustering samples
################################################################################
source('analysis/09-colors.R')

library(pvclust)
library(dendextend)

set.seed(1234)

tmat <- as.matrix(assay(tform.sel))
samples.pv <- pvclust::pvclust(tmat,
                               method.dist = 'correlation',
                               use.cor = "pairwise.complete.obs",
                               method.hclust = 'average', 
                               nboot=1000
)

if(makepdf) pdf('analysis/figures/bc.sample_clustering.pdf', 
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
    plot(main="Best Counts - Cluster dendrogram with AU/BP values")

samples.pv %>% text(cex=0.7) 

colored_dots(data.frame(
    "paired"=layout_colors[as.character(colData(tform.sel)[dorder,'layout'])],
    "center"=center_colors[as.character(colData(tform.sel)[dorder,'sequencing_center'])],
    "cancer"=karyotype_colors[as.character(colData(tform.sel)[dorder,'karyotype'])],
    celltype=dcols
), y_shift=-0.14)

if(makepdf) dev.off()

################################################################################
# Best counts - Differentially expressed genes
################################################################################
# Loads dend.short and celltypes_order
load('analysis/Rdata/pvclust_celltypes.Rdata')

# Adjusted p-value (FDR) cutoff
padj_cutoff <- 0.1

# log2FoldChange (LFC) cutoff
lfc_cutoff <- 1.0

################################################################################
# DESeq2 testing results
################################################################################
ntypes <- length(celltypes_order)

pwresult <- lapply(1:ntypes, function(i) {
    ret <- lapply(1:ntypes, function(j){
        if(i == j) return(NULL)
        results(dds.sel, contrast=c('celltype', celltypes_order[j], celltypes_order[i]))
    })
    names(ret) <- celltypes_order
    ret
})
names(pwresult) <- celltypes_order

degcounts <- lapply(pwresult, function(reslist) {
    sapply(reslist, function(v) {
        if(is.null(v)) return(0)
        subset(v, padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) %>% nrow
    })
}) %>% bind_cols() %>% data.frame
row.names(degcounts) <- colnames(degcounts) <- names(pwresult)

if(makepdf) pdf('analysis/figures/bcounts.celltype_deg.pdf', 
    paper='USr', width=11, height=8.5, onefile = FALSE)

hmbreaks.deg <- c(0, seq(5, 1500, 5), max(degcounts, 1501))
hmcolor.deg <- c('#FFFFFF', colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(hmbreaks.deg)-2))

pheatmap::pheatmap(degcounts,
                   color=hmcolor.deg,
                   breaks=hmbreaks.deg,
                   legend_breaks=c(0,500,1000,1500),
                   display_numbers=TRUE, number_format='%d', 
                   cluster_rows=FALSE, cluster_cols=FALSE,
                   cellwidth=20, cellheight=20)

if(makepdf) dev.off()

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
    geom_point(size=3) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
    scale_color_manual(values=celltype_colors) +
    ggtitle('Best Counts - HERV expression PCA')

if(makepdf) pdf('analysis/figures/bcounts.sample_pca_1v2.pdf', paper='USr', width=10, height=8)
p12
if(makepdf) dev.off()

p13 <- pca.df  %>%
    ggplot(aes_string(x = "PC1", y = "PC3", color = "celltype", shape='sequencing')) + 
    geom_point(size=2) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
    scale_color_manual(values=celltype_colors) +
    ggtitle('Best Counts - HERV expression PCA')

if(makepdf) pdf('analysis/figures/bcounts.sample_pca_1v3.pdf', paper='USr', width=10, height=8)
p13
if(makepdf) dev.off()

################################################################################
# Save data
################################################################################
best_counts <- list(
    'dds.sel' = dds.sel,
    'tform.sel' = tform.sel,
    'tform.sel.0' = tform.sel.0,
    'sfilt' = sfilt,
    'gfilt' = gfilt,
    'samples.pv' = samples.pv,
    'pwresult' = pwresult
)

save(best_counts, file='analysis/Rdata/best_counts.Rdata')
