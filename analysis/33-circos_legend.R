#! /usr/bin/env Rscript
library(tidyverse)
source('analysis/09-colors.R')

circos_dir <- 'analysis/circos'
ok_chroms <- c(sprintf('chr%d', 1:22), 'chrX', 'chrY')

################################################################################
# HERV counts by windows
################################################################################
winsize <- 10000000

#--- Get chromosome lengths
chromlengths <- 
    data.frame(rawtext=readLines(file.path('refs', 'hg38_header.txt')), 
               stringsAsFactors = F) %>%
    tidyr::separate(rawtext, c('sq', 'f1', 'f2'), sep='\t', remove=FALSE) %>% 
    tidyr::separate(f1, c('sn', 'chrom'), sep=':') %>% 
    tidyr::separate(f2, c('ln', 'len'), sep=':') %>%
    mutate(len=as.numeric(len)) %>%
    select(chrom, len) %>%
    dplyr::filter(chrom %in% ok_chroms)

clens <- chromlengths$len
names(clens) <- chromlengths$chrom

x4 <- seq(1, clens[4], length.out=20)
tmp4 <- data.frame(chrom='hs4', start=x4[1:(length(x4)-1)], end=x4[2:length(x4)], stringsAsFactors = F)
width <- (tmp4$end - tmp4$start)[1]

x3 <- seq(clens[3], 1, -width)
tmp3 <- data.frame(chrom='hs3', start=rev(x3[2:length(x3)]), end=rev(x3[1:(length(x3)-1)]),
                   stringsAsFactors = F)

x5 <- seq(1, clens[5], width)
tmp5 <- data.frame(chrom='hs5', start=x5[1:(length(x5)-1)], end=x5[2:length(x5)],
                   stringsAsFactors = F)

legend <- rbind(tmp3, tmp4, tmp5)
legend <- legend[3:54,]
legend$colorname <- c('heatcolor-min', sprintf('heatcolor-50-seq-%d', 1:50), 'heatcolor-max')
legend$value <- c(0, seq(0.1, 9.9, by=0.2), 10)
legend$end <- legend$end - 1

outfile <- file.path(circos_dir, 'data', 'heatmap_legend.txt')
legend %>% 
    dplyr::mutate(
        start=as.integer(start),
        end=as.integer(end),
        opts=sprintf('color=%s', colorname)
    ) %>%
    select(chrom, start, end, value, opts) %>%
    write.table(file=outfile, sep='\t', row.names=F, col.names=F, quote=F)

outfile <- file.path(circos_dir, 'data', 'heatmap_legend_text.txt')
legend[c(1,26,52),] %>%
    dplyr::mutate(
        start=as.integer(start),
        end=as.integer(end),
        value=round(value)
    ) %>%
    select(chrom, start, end, value) %>%
    write.table(file=outfile, sep='\t', row.names=F, col.names=F, quote=F)
