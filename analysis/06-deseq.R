#! /usr/bin/env Rscript
load('analysis/Rdata/loaded.Rdata')

library(tidyverse)
library(DESeq2)


do.parallel <- FALSE
if(do.parallel) {
    library(BiocParallel)
    register(MulticoreParam(8))
}

################################################################################
# All samples
################################################################################
dds.all <- DESeqDataSetFromMatrix(counts.telescope, samples, ~celltype)
sizeFactors(dds.all) <- calculateSizeFactors(metrics$mapped_frags.bt2)
dds.all <- DESeq(dds.all, parallel=do.parallel)
tform.all <- varianceStabilizingTransformation(dds.all, blind=FALSE)
tform.all.0 <- varianceStabilizingTransformation(dds.all, blind=TRUE)

save(dds.all, tform.all, tform.all.0, file='analysis/Rdata/deseq.all.Rdata')

################################################################################
# PolyA extract only
################################################################################
# Criteria for samples
extracts <- c("longPolyA")
sfilt <- samples$rnaextract %in% extracts

# Criteria for genes
mincount <- 5
minsamp  <- 2
gfilt <- rowSums(counts.telescope[,sfilt] >= mincount) >= minsamp
gfilt <- gfilt & (annot.herv$chrom != 'chrY')

dds.sel <- DESeqDataSetFromMatrix(counts.telescope[gfilt, sfilt], samples[sfilt, ], ~celltype)
sizeFactors(dds.sel) <- calculateSizeFactors(metrics[sfilt,]$mapped_frags.bt2)
dds.sel <- DESeq(dds.sel, parallel=do.parallel)
tform.sel <- varianceStabilizingTransformation(dds.sel, blind=FALSE)
tform.sel.0 <- varianceStabilizingTransformation(dds.sel, blind=TRUE)

save(dds.sel, tform.sel, tform.sel.0, sfilt, gfilt,
     file='analysis/Rdata/deseq.sel.Rdata')
