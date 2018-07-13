#! /usr/bin/env Rscript
if(!exists('counts.telescope')) source('analysis/03-load_counts.R')
if(!exists('metrics')) source('analysis/04-load_metrics.R')

library(tidyverse)

################################################################################
# Counts Per Million
################################################################################
stopifnot(all(colnames(counts.telescope) == rownames(metrics)))

# CPM for telescope, best, and TETranscripts - bowtie2 mapping
cpm.telescope <- t(t(counts.telescope) / (metrics$mapped_frags.bt2 / 1e6))
cpm.best <- t(t(counts.best) / (metrics$mapped_frags.bt2 / 1e6))
cpm.fam <- t(t(counts.fam) / (metrics$mapped_frags.bt2 / 1e6))

# CPM for unique - only uniquely mapping fragments
cpm.unique <- t(t(counts.unique) / (metrics$unique.ts / 1e6))

# CPM for RepEnrich - bowtie1 mapped frags
cpm.repenrich <- t(t(counts.repenrich) / (metrics$mapped_frags.bt1 / 1e6))

# CPM for TETranscripts - bowtie2 mapped frags
cpm.tetx <- t(t(counts.tetx) / (metrics$mapped_frags.bt2 / 1e6))
