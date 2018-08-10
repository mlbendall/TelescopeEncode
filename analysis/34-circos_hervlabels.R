#! /usr/bin/env Rscript
source('analysis/07-expressed_cpm.R')

library(tidyverse)

circos_dir <- 'analysis/circos'
ok_chroms <- c(sprintf('chr%d', 1:22), 'chrX', 'chrY')

################################################################################
# 
################################################################################
annot_exp <- annot.herv %>% 
    cbind(exp.any=expressed.celltype.any, exp.all=expressed.celltype.all) %>%
    cbind(expressed.celltype)

outfile <- file.path(circos_dir, 'data', 'herv.text.txt')
annot_exp %>%
    dplyr::filter(chrom %in% ok_chroms & exp.any) %>%
    dplyr::mutate(
        chrom=gsub('^chr', 'hs', chrom),
        value=locus
    ) %>%
    select(chrom, start, end, value) %>%
    write.table(file=outfile, sep='\t', row.names=F, col.names=F, quote=F)
