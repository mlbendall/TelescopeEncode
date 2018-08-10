#! /usr/bin/env Rscript
load('analysis/Rdata/loaded.Rdata')

library(tidyverse)

################################################################################
# Select samples to consider for "expressed"
################################################################################
# extracts <- c("longPolyA")
# 
# celltypes <- levels(samples$celltype) # All celltypes
# 
# # Determine which celltypes have replicates
# samp_cts <- table(samples[samples$rnaextract %in% extracts,]$celltype)
# ctypes_rep <- names(samp_cts[samp_cts > 1])
# rm(samp_cts)

################################################################################
# Counts Per Million
################################################################################
mincpm <- 0.5

# Expression threshold for each sample
expressed.sample <- cpm.telescope > mincpm

# Considering only polyA samples, locus expressed in at least `minsamp` samples
extracts <- c('longPolyA')
minsamp  <- 2
expressed.core <- rowSums(expressed.sample[ , samples$rnaextract %in% extracts]) >= minsamp

# Locus expressed per cell type in majority of replicates
expressed.celltype <- lapply(celltypes.all, function(ctype) {
    sfilt <- (samples$rnaextract %in% extracts) & (samples$celltype == ctype)
    tmp.minsamp <- sum(sfilt) / 2
    (rowSums(expressed.sample[, sfilt, drop=FALSE]) > tmp.minsamp) %>%
        data.frame
}) %>% bind_cols
colnames(expressed.celltype) <- celltypes.all
rownames(expressed.celltype) <- row.names(cpm.telescope)

colSums(expressed.celltype)

expressed.celltype.any <- rowSums(expressed.celltype[,celltypes.rep]) > 0
expressed.celltype.all <- rowSums(expressed.celltype[,celltypes.rep]) == length(celltypes.rep)


