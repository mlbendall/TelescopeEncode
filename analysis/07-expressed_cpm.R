#! /usr/bin/env Rscript
load('analysis/Rdata/loaded.Rdata')

library(tidyverse)

################################################################################
# Select samples to consider for "expressed"
################################################################################
extracts <- c("longPolyA")

celltypes <- levels(samples$celltype) # All celltypes

samp_cts <- table(samples[samples$rnaextract %in% extracts,]$celltype)
ctypes_rep <- names(samp_cts[samp_cts > 1])
rm(samp_cts)

################################################################################
# Counts Per Million
################################################################################
mincpm <- 0.5
minsamp  <- 2

expressed.core <- rowSums(cpm.telescope[ , samples$rnaextract %in% extracts] > mincpm) >= minsamp

expressed.celltype <- lapply(celltypes, function(ctype) {
    sfilt <- (samples$rnaextract %in% extracts) & (samples$celltype == ctype)
    ret <- data.frame(rowSums(cpm.telescope[, sfilt, drop=FALSE] > mincpm) == sum(sfilt))
    names(ret) <- ctype
    ret
}) %>% bind_cols
rownames(expressed.celltype) <- row.names(cpm.telescope)

expressed.celltype.any <- rowSums(expressed.celltype[,ctypes_rep]) > 0
expressed.celltype.all <- rowSums(expressed.celltype[,ctypes_rep]) == length(ctypes_rep)
