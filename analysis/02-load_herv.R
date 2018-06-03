#! /usr/bin/env Rscript

library(tidyverse)

## Load HERV annotation

# Load mapping between HERV families and nicknames. This is used for coloring
# and visualization.
herv_fam <- read.table('analysis/herv_families.tsv', 
                       sep='\t', header=T, stringsAsFactors = F)

# Load TSV file with one row for each HERV in the annotation. Telescope reports
# are put into the same order as this TSV, and nicknames are added to the table.
annot.herv <- read.table('refs/HERV_rmsk.hg38.v2.tsv',
                         sep='\t', header=T, stringsAsFactors = F) %>%
    dplyr::mutate(length=end - start + 1) %>%
    tidyr::separate(locus, c("family"), sep='_', extra='drop', remove=F) %>%
    dplyr::left_join(herv_fam, by='family') %>%
    dplyr::mutate(
        family=factor(family, levels=herv_fam$family),
        group=factor(group, levels=unique(herv_fam$group)),
        letter=factor(letter, levels=unique(herv_fam$letter))
    ) %>%
    dplyr::select(locus, chrom, start, end, length, family, group, letter, category)

# Final Counts from Telescope
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

# Unique Counts calculated by Telescope
counts.unique <- lapply(samples$sample,
                     function(s){
                         tmp <- read.table(file.path('samples', s, 'inform-telescope_report.tsv'),
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
row.names(counts.unique) <- annot.herv$locus

# Best Counts calculated by Telescope
counts.best <- lapply(samples$sample,
                        function(s){
                            tmp <- read.table(file.path('samples', s, 'inform-telescope_report.tsv'),
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
row.names(counts.best) <- annot.herv$locus
