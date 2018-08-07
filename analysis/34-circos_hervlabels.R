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

outfile <- file.path(circos_dir, 'data', 'herv.textB.txt')
    
annot_exp %>%
    dplyr::filter(chrom == 'chr19' & start > 53000000 & end < 59000000 & exp.any) %>%
    dplyr::mutate(
        chrom=gsub('^chr', 'hs', chrom),
        value=locus
    ) %>%
    select(chrom, start, end, value) %>%
    write.table(file=outfile, sep='\t', row.names=F, col.names=F, quote=F)


    # band hs19 q13.41 q13.41 50900000 53100000 gpos25
    # band hs19 q13.42 q13.42 53100000 55800000 gneg
    # band hs19 q13.43 q13.43 55800000 58617616 gpos25

outfile <- file.path(circos_dir, 'data', 'herv.textC.txt')
annot_exp %>%
    dplyr::filter(chrom == 'chr6' & end < 60000000 & exp.any) %>%
    dplyr::mutate(
        chrom=gsub('^chr', 'hs', chrom),
        value=locus
    ) %>%
    select(chrom, start, end, value) %>%
    write.table(file=outfile, sep='\t', row.names=F, col.names=F, quote=F)
