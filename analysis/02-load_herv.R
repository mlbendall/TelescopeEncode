#! /usr/bin/env Rscript

library(tidyverse)
# library(tximport)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(12))

## Load samples
if(file.exists('load_sampledata.Rdata')) {
    load('load_sampledata.Rdata')
} else {
    source('load_sampledata.R')
}
rm(list=ls(pattern='samp\\.'))
rm(celltypes, lineages, stages, studies, study_abbr)



## Load HERV annotation
herv_fam <- read.table('analysis/herv_families.tsv', sep='\t', header=T, stringsAsFactors = F)
annot.herv <- read.table('refs/HERV_rmsk.hg38.v2.tsv.gz',
                         sep='\t', header=T, stringsAsFactors = F) %>%
    dplyr::mutate(length=end - start + 1) %>%
    tidyr::separate(locus, c("family"), sep='_', extra='drop', remove=F) %>%
    dplyr::left_join(herv_fam, by='family') %>%
    dplyr::mutate(
        family=factor(family, levels=herv_fam$family),
        group=factor(group, levels=unique(herv_fam$group))
    ) %>%
    dplyr::select(locus, chrom, start, end, length, family, group, category)

### Remove sample(s) that failed for telescope
samples[!(row.names(samples) %in% c('GSM958733')),] %>%
    droplevels -> samples

### Final Counts from Telescope
counts.herv <- lapply(samples$sample_id,
                      function(s){
                          tmp <- read.table(file.path('results', paste(s, 'inform-telescope_report.tsv', sep='.')),
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
row.names(counts.herv) <- annot.herv$locus

### Unique Counts from Telescope
counts.unique <- lapply(samples$sample_id,
                      function(s){
                          tmp <- read.table(file.path('results', paste(s, 'inform-telescope_report.tsv', sep='.')),
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

## Create DESeq object
dds.rx <- DESeqDataSetFromMatrix(counts.herv, samples, ~1)
dds.rx <- DESeq(dds.rx, parallel=F)

## Transform DESeq counts
tform.rx.vsd <- varianceStabilizingTransformation(dds.rx, blind=FALSE)
tform.rx.norm <- normTransform(dds.rx)

