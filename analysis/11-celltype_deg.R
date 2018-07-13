#! /usr/bin/env Rscript

load('analysis/Rdata/deseq.sel.Rdata')
load('analysis/Rdata/pvclust_celltypes.Rdata')
source('analysis/09-colors.R')

library(tidyverse)
library(DESeq2)

makepdf <- TRUE

################################################################################
# Parameters
################################################################################
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

save(pwresult, file='analysis/Rdata/pwresult_deg.Rdata')

degcounts <- lapply(pwresult, function(reslist) {
    sapply(reslist, function(v) {
        if(is.null(v)) return(0)
        subset(v, padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) %>% nrow
    })
}) %>% bind_cols() %>% data.frame
row.names(degcounts) <- colnames(degcounts) <- names(pwresult)

if(makepdf) pdf('analysis/figures/celltype_deg.pdf', 
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
