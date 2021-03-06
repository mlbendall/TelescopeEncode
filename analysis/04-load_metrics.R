#! /usr/bin/env Rscript
if(!exists('samples')) source('analysis/01-load_sampledata.R')

library(tidyverse)

# Function for calculating size factors
calculateSizeFactors <- function(mapped_frags) {
    geomean <- expm1(mean(log1p(mapped_frags)))
    mapped_frags / geomean
}


################################################################################
# Load metrics data
################################################################################

#--- Load bowtie2 alignment metrics
met.bt2 <- lapply(samples$sample,
               function(s){
                   tmp <- read.table(file.path('samples', s, 'bt2_multi.summary.txt'),
                              header=T, stringsAsFactors = F) %>%
                       filter(run == "total")
                   if( 'unpaired' %in% names(tmp))
                       tmp <- tmp %>%
                           mutate(concord_0=0, concord_1=0, concord_M=0, discord_1=0)
                   tmp %>%
                       mutate(sample=s) %>%
                       select(sample, reads, concord_0, concord_1, concord_M, 
                              discord_1, aln_0, aln_1, aln_M, arate)
               }) %>% 
    bind_rows() %>% data.frame(stringsAsFactors=F)

names(met.bt2)[2:10] <- paste(names(met.bt2)[2:10], 'bt2', sep='.')

#--- Load bowtie1 alignment metrics
met.bt1 <- lapply(samples$sample,
                  function(s){
                      read.table(file.path('samples', s, 'bt_repenrich.summary.txt'),
                                        header=T, stringsAsFactors = F) %>%
                          filter(run == "total") %>%
                          mutate(sample=s) %>%
                          select(sample, reads, aln_0, aln_1, aln_M, pct_0, pct_1=pct, pct_M)

                  }) %>%
    bind_rows() %>% data.frame(stringsAsFactors=F)

names(met.bt1)[2:8] <- paste(names(met.bt1)[2:8], 'bt1', sep='.')

#--- Load telescope metrics from comment
metrics.list <- lapply(samples$sample, 
                       function(s){
                           f <- file.path('samples', s, 'inform-telescope_report.tsv')
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
met.ts <- lapply(metrics.list, 
                 function(m) {
                     ret <- sapply(mn, function(x) m[x])
                     names(ret) <- gsub('\\..*', '', names(ret))
                     ret
                 }) %>% do.call(rbind, .) %>% data.frame
names(met.ts) <- paste(names(met.ts), 'ts', sep='.')
met.ts$sample <- samples$sample

#--- Merge all metrics
metrics <- dplyr::inner_join(met.bt1, met.bt2, by='sample') %>%
    dplyr::inner_join( met.ts, by='sample')

row.names(metrics) <- metrics$sample

#--- Add library size (total frags mapped) columns
metrics$mapped_frags.bt2 <- metrics$pair_mapped.ts + metrics$pair_mixed.ts + metrics$single_mapped.ts
metrics$mapped_frags.bt1 <- metrics$reads.bt1 - metrics$aln_0.bt1

rm(met.ts, met.bt2, met.bt1, mn, metrics.list)
