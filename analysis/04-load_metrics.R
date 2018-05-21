#! /usr/bin/env Rscript

library(tidyverse)

## Load samples
if(file.exists('load_sampledata.Rdata')) {
    load('load_sampledata.Rdata')
} else {
    source('load_sampledata.R')
}
rm(list=ls(pattern='samp\\.'))
rm(celltypes, lineages, stages, studies, study_abbr)


# Load alignment metrics
met.aln <- lapply(samples$sample,
               function(s){
                   lines <- readLines(file.path('genequant', s, 'bt2_multi.summary.txt'))
                   c(s,
                     as.numeric(gsub('^(\\d+) .*', '\\1',  lines[1])),
                     as.numeric(gsub('^(\\d+\\.\\d+)% .*', '\\1',  lines[length(lines)]))
                   )
               }) %>% do.call(rbind, .) %>% data.frame(stringsAsFactors=F)

names(met.aln) <- c('sample', 'total_reads', 'alnrate')
met.aln$total_reads <- as.integer(met.aln$total_reads)
met.aln$alnrate <- as.double(met.aln$alnrate) * 1e-2
met.aln <- met.aln %>% mutate(naln=floor(total_reads*alnrate))
row.names(met.aln) <- met.aln$sample


# Load telescope metrics from comment
metrics.list <- lapply(samples$sample, 
                       function(s){
                           f <- file.path('genequant', s, 'inform-telescope_report.tsv')
                           if(file.exists(f)) {
                               h <- readLines(f, 1) %>% strsplit(., '\t') %>% unlist
                               rstr <- sapply(strsplit(h[-c(1,2)], ':'), function(t) as.numeric(unlist(t[2][1])))
                               names(rstr) <- sapply(strsplit(h[-c(1,2)], ':'), function(t) t[1])
                           } else{
                               rstr <- c(NA)
                           }
                           rstr
                       }
)
mn <- unique(do.call(c, lapply(metrics.list, names)))
met.ts <- lapply(metrics.list, function(m) {
    ret <- sapply(mn, function(x) m[x])
    names(ret) <- gsub('\\..*', '', names(ret))
    ret
}) %>% do.call(rbind, .) %>% data.frame
row.names(met.ts) <- samples$sample

## Save Rdata
save(samples, met.aln, met.ts, file='load_metrics.Rdata')
