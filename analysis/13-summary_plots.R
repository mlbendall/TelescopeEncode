#! /usr/bin/env Rscript
source('analysis/07-expressed_cpm.R')

library(tidyverse)
library(cowplot)

source('analysis/09-colors.R')

makepdf <- TRUE

################################################################################
# Proportion of fragments mapping to HERV
################################################################################
stopifnot(all(colnames(cpm.telescope) == row.names(samples)))

tmp <- data.frame(cpm=colSums(cpm.telescope)) %>%
    tibble::rownames_to_column('sample') %>%
    dplyr::left_join(samples, by='sample') %>%
    dplyr::filter(rnaextract %in% extracts & celltype %in% celltypes.rep) %>%
    dplyr::mutate(celltype=factor(celltype, levels=celltypes.rep)) %>%
    dplyr::mutate(sequencing=paste(sequencing_center, layout)) %>%
    dplyr::mutate(sequencing=factor(sequencing, levels=c('CSHL PAIRED', 'CALTECH SINGLE', 'CALTECH PAIRED'))) %>%
    dplyr::mutate(prop=(cpm/1e6)) %>%
    dplyr::select(sample, cpm, prop, celltype, sequencing)

p <- tmp %>%
    ggplot(aes(x=celltype, y=prop)) +
    geom_boxplot(aes(fill=celltype), outlier.alpha=0, color='#333333') +
    geom_point(size=2, color='#333333') +
    scale_fill_manual(values=celltype_colors) +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') + ylab('Proportion mapping to HERV')

if(makepdf) pdf('analysis/figures/herv_proportion.pdf', paper='USr', width=10, height=5)
p
if(makepdf) dev.off()

################################################################################
# Number of HERV elements expressed per celltype
################################################################################
sel <- expressed.celltype[,celltypes.rep]
specific <- lapply(celltypes.rep, function(ctype) {
    sel[,ctype] & rowSums(sel) == 1
}) %>% do.call(cbind, .) %>% data.frame()
colnames(specific) <- celltypes.rep
row.names(specific) <- row.names(sel)

p <- data.frame(expressed=colSums(sel), spec=colSums(specific)) %>%
    tibble::rownames_to_column('celltype') %>%
    dplyr::mutate(shared=expressed-spec) %>%
    dplyr::select(-expressed) %>%
    tidyr::gather('key', 'count', 2:3) %>%
    dplyr::mutate(
        celltype=factor(celltype, levels=celltypes.rep),
        key=factor(key, levels=c('spec', 'shared'))
    ) %>%
    ggplot(aes(x=celltype, y=count, fill=celltype, alpha=key)) +
    geom_col(color='#333333') +
    scale_fill_manual(values=celltype_colors) +
    scale_alpha_discrete(name='', labels=c('Cell-type specific', 'Shared'), range=c(1.0, 0.60)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') + ylab('# HERV elements expressed')

if(makepdf) pdf('analysis/figures/herv_elements.pdf', paper='USr', width=10, height=5)
p
if(makepdf) dev.off()

