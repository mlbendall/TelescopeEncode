#! /usr/bin/env Rscript
library(tidyverse)

if(!exists('metrics')) source('analysis/01-load_sampledata.R')

################################################################################
# Load metrics data
################################################################################

library(tidyverse)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(12))



# Calculate size factors from library size
mapped_frags <- metrics$pair_mapped + metrics$pair_mixed + metrics$single_mapped
geomean <- expm1(mean(log1p(mapped_frags)))
sizefactors.gm <- mapped_frags / geomean

# Calculate CPM
cpm.telescope <- t(t(counts.telescope) / (mapped_frags / 1e6))

# Create DESeq object
dds <- DESeqDataSetFromMatrix(counts.telescope, samples, ~celltype)
sizeFactors(dds) <- sizefactors.gm
dds <- DESeq(dds, parallel=F)

# Transform DESeq counts
tform <- varianceStabilizingTransformation(dds, blind=FALSE)

