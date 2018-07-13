#! /usr/bin/env Rscript
if(!exists('samples')) source('analysis/01-load_sampledata.R')
if(!exists('annot.herv')) source('analysis/02-load_herv.R')

library(tidyverse)

################################################################################
# Load count data
################################################################################

#--- Final Counts from Telescope
counts.telescope <- lapply(samples$sample,
                           function(s){
                               tmp <- read.table(
                                   file.path('samples', s, 'inform-telescope_report.tsv'),
                                   sep='\t', header=T, stringsAsFactors=F)
                               ret <- data.frame(transcript=annot.herv$locus, stringsAsFactors=F) %>%
                                   left_join(tmp, by='transcript') %>%
                                   mutate(
                                       gene_id = transcript,
                                       count = final_count
                                   ) %>%
                                   select(gene_id, count)
                               ret[is.na(ret)] <- 0
                               stopifnot(all(ret$gene_id == annot.herv$locus))
                               ret$gene_id <- NULL
                               names(ret) <- c(s)
                               ret
                           }) %>%
    bind_cols
row.names(counts.telescope) <- annot.herv$locus

#--- Telescope summarized at family level
counts.fam <- counts.telescope %>%
    dplyr::mutate(family=annot.herv$family) %>%
    dplyr::group_by(family) %>%
    dplyr::summarize_if(is.numeric, sum) %>%
    data.frame

row.names(counts.fam) <- counts.fam$family
counts.fam$family <- NULL


#--- Unique Counts calculated by Telescope
counts.unique <- lapply(samples$sample,
                        function(s){
                            tmp <- read.table(file.path('samples', s, 'inform-telescope_report.tsv'),
                                              sep='\t', header=T, stringsAsFactors=F)
                            ret <- data.frame(transcript=annot.herv$locus, stringsAsFactors=F) %>%
                                left_join(tmp, by='transcript') %>%
                                mutate(
                                    gene_id = transcript,
                                    count = unique_count
                                ) %>%
                                select(gene_id, count)
                            ret[is.na(ret)] <- 0
                            stopifnot(all(ret$gene_id == annot.herv$locus))
                            ret$gene_id <- NULL
                            names(ret) <- c(s)
                            ret
                        }) %>%
    bind_cols
row.names(counts.unique) <- annot.herv$locus

#--- Best Counts calculated by Telescope
counts.best <- lapply(samples$sample,
                      function(s){
                          tmp <- read.table(file.path('samples', s, 'inform-telescope_report.tsv'),
                                            sep='\t', header=T, stringsAsFactors=F)
                          ret <- data.frame(transcript=annot.herv$locus, stringsAsFactors=F) %>%
                              left_join(tmp, by='transcript') %>%
                              mutate(
                                  gene_id = transcript,
                                  count = init_best_random
                              ) %>%
                              select(gene_id, count)
                          ret[is.na(ret)] <- 0
                          stopifnot(all(ret$gene_id == annot.herv$locus))
                          ret$gene_id <- NULL
                          names(ret) <- c(s)
                          ret
                      }) %>%
    bind_cols
row.names(counts.best) <- annot.herv$locus

#--- Fraction counts from RepEnrich
REfams <- read.table(file.path('samples', samples$sample[1], 'RepEnrich_fraction_counts.txt'),
                     sep='\t', header=F, stringsAsFactors = F) %>%
           select(family=V1)

counts.repenrich <- lapply(samples$sample,
       function(s) {
           read.table(file.path('samples', s, 'RepEnrich_fraction_counts.txt'),
                      sep='\t', header=F, stringsAsFactors = F) %>%
           select(V4)
       }) %>% bind_cols()

names(counts.repenrich) <- samples$sample
row.names(counts.repenrich) <- REfams$family
rm(REfams)


#--- Counts from TETranscripts
last_gene <- 'gSpikein_phiX174'

# Get list of families
tmp <- read.table(file.path('samples', samples$sample[1], 'tetx_count.cntTable'),
                  sep='\t', header=T, stringsAsFactors=F)
te <- tmp[(which(tmp$gene.TE == last_gene) + 1):nrow(tmp), ]
REfams <- te %>%
    tidyr::separate(gene.TE, c('family', 'superfamily', 'class'), sep=':') %>%
    filter(class == 'LTR') %>%
    select(family)

rm(tmp, te)

counts.tetx <- lapply(samples$sample,
    function(s) {
        tmp <- read.table(file.path('samples', s, 'tetx_count.cntTable'),
                          sep='\t', header=T, stringsAsFactors=F)
        te <- tmp[(which(tmp$gene.TE == last_gene) + 1):nrow(tmp), ]
        te <- te %>%
            tidyr::separate(gene.TE, c('family', 'superfamily', 'class'), sep=':') %>%
            filter(class == 'LTR')
        te <- dplyr::left_join(REfams, te, by='family')
        names(te)[4] <- s
        te %>% dplyr::select(s)
    }) %>% bind_cols

row.names(counts.tetx) <- REfams$family
rm(REfams, last_gene)
