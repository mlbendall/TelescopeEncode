#! /usr/bin/env Rscript

source('analysis/01-load_sampledata.R')
source('analysis/02-load_herv.R')
source('analysis/03-load_counts.R')
source('analysis/04-load_metrics.R')
source('analysis/05-calculate_cpm.R')

save.image('analysis/Rdata/loaded.Rdata')
