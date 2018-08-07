#! /usr/bin/env Rscript
source('analysis/07-expressed_cpm.R')

library(tidyverse)

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

lapply(1:nrow(chromlengths), function(r){
    ret <- data.frame(chrom=chromlengths[r,1], 
                      winstart=seq(0, chromlengths[r,2], winsize),
                      stringsAsFactors = F)
    ret$winend <- pmin(chromlengths[r,2], ret$winstart + winsize)
    ret
} ) %>% dplyr::bind_rows() -> gwindows

tmp <- annot.herv %>%
    cbind(exp.any=expressed.celltype.any, exp.all=expressed.celltype.all) %>%
    mutate(midpoint=(start + end) / 2)

tmpwin <- lapply(1:nrow(gwindows), function(r) {
    inwin <- tmp %>%
        dplyr::filter(chrom == gwindows[r,]$chrom) %>% 
        dplyr::filter(midpoint >= gwindows[r,]$winstart) %>%
        dplyr::filter(midpoint < gwindows[r,]$winend)
    c(count=nrow(inwin),
      expressed=sum(inwin$exp.any),
      table(inwin$family),
      table(inwin$letter)
    )
}) %>% do.call(rbind, .) %>% data.frame

gwindows <- cbind(gwindows, tmpwin)

#--- All HERV loci
gwindows %>%
    mutate(
        chr=gsub('^chr', 'hs', chrom),
        winstart=as.integer(winstart),
        winend=as.integer(winend),
    ) %>%
    select(chr, winstart, winend, count) %>%
    write.table(
        file=file.path(circos_dir, 'data', 'herv_count.histogram.txt'),
        sep='\t', row.names=F, col.names=F, quote=F
    )

#--- Expressed HERV loci
gwindows %>%
    mutate(
        chr=gsub('^chr', 'hs', chrom),
        winstart=as.integer(winstart),
        winend=as.integer(winend),
    ) %>%
    select(chr, winstart, winend, expressed) %>%
    write.table(
        file=file.path(circos_dir, 'data', 'herv_expressed.histogram.txt'),
        sep='\t', row.names=F, col.names=F, quote=F
    )

#--- Silent HERV loci
gwindows %>%
    mutate(
        chr=gsub('^chr', 'hs', chrom),
        winstart=as.integer(winstart),
        winend=as.integer(winend),
        silent=count - expressed
    ) %>%
    select(chr, winstart, winend, silent) %>%
    write.table(
        file=file.path(circos_dir, 'data', 'herv_silent.histogram.txt'),
        sep='\t', row.names=F, col.names=F, quote=F
    )

#--- Generate windowhist.conf
r1_start <- 0.90
hist_start <- 0.995

outfile <- file.path(circos_dir, 'windowhist.conf')
cat(
    '<plot>',
    'type=histogram',
    'file=data/herv_count.histogram.txt',
    sprintf('r1 = %0.3fr', hist_start),
    sprintf('r0 = %0.3fr', r1_start + 0.005),
    'min = 0',
    'max = 200',
    'stroke_type = bin',
    'thickness   = 2',
    'color       = vdgrey',
    'fill_color  = red',
    'extend_bin  = no',
    '<axes>', '<axis>',
    'spacing   = 0.25r',
    'color     = lgrey',
    'thickness = 2',
    '</axis>','</axes>',
    '</plot>', 
    '',
    '<plot>',
    'type=histogram',
    'file=data/herv_silent.histogram.txt',
    sprintf('r1 = %0.3fr', hist_start),
    sprintf('r0 = %0.3fr', r1_start + 0.005),
    'min = 0',
    'max = 200',
    'stroke_type = bin',
    'thickness   = 2',
    'color       = vdgrey',
    'fill_color  = vvlgrey',
    'extend_bin  = no',
    '</plot>',
    sep='\n',
    file=outfile, append=FALSE
)
