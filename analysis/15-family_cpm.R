#! /usr/bin/env Rscript
source('analysis/07-expressed_cpm.R')

library(tidyverse)
library(cowplot)

source('analysis/09-colors.R')

makepdf <- TRUE

################################################################################
# Selected families
################################################################################
selfams <- c('HERVH', 'HERVE', 'HML2', 'HML6', 
             'HERV9', 'HERVL', 'HARLEQUIN', 'HERV3')



################################################################################
# Expressed HERV elements (consensus), celltype X family
################################################################################

# Consensus expression (logical) for each locus, selected celltypes
expressed.celltype.sel <- expressed.celltype[,celltypes.rep]

# Locus count per family, per celltype
loccount.celltype <- expressed.celltype.sel %>%
    dplyr::mutate(family=annot.herv$family) %>%
    dplyr::group_by(family) %>%
    dplyr::summarize_if(is.logical, sum) %>%
    tidyr::gather('celltype', 'count', -c('family')) %>%
    dplyr::mutate(celltype=factor(celltype, levels=celltypes.rep)) %>%
    dplyr::select(family,count,celltype)

p.all <- loccount.celltype %>%
    ggplot(aes(celltype, count, label=count)) +
    geom_col(aes(fill=celltype)) +
    geom_text(size = 3, position = position_stack(vjust = 1.07), color='#010101') +   
    scale_fill_manual(values=celltype_colors, guide=FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~family, ncol=1, scales='free_y', strip.position='right') +
    xlab('') + ylab('# HERV elements expressed')

if(makepdf) pdf('analysis/figures/loccount_ALL.celltypeXfamily.pdf', width=4, height=120)
p.all
if(makepdf) dev.off()

p <- loccount.celltype %>%
    dplyr::filter(family %in% selfams) %>%
    dplyr::mutate(family=factor(family, levels=selfams)) %>%
    ggplot(aes(celltype, count, label=count)) +
    geom_col(aes(fill=celltype)) +
    geom_text(size = 3, position = position_stack(vjust = 1.07), color='#010101') +   
    scale_fill_manual(values=celltype_colors, guide=FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~family, ncol=1, scales='free_y', strip.position='right') +
    xlab('') + ylab('# HERV elements expressed')

if(makepdf) pdf('analysis/figures/loccount.celltypeXfamily.pdf', paper='letter', width=4, height=11)
p
if(makepdf) dev.off()

################################################################################
# Same as above but shade specific vs. shared
################################################################################
# Celltype-specific expression (logical) for each locus
specific <- lapply(celltypes.rep, function(ctype) {
    expressed.celltype.sel[,ctype] & rowSums(expressed.celltype.sel) == 1
}) %>% do.call(cbind, .) %>% data.frame()
colnames(specific) <- celltypes.rep

# Celltype-specific locus count per family
loccount.celltype.specific <- specific %>% 
    dplyr::mutate(family=annot.herv$family) %>%
    dplyr::group_by(family) %>%
    dplyr::summarize_if(is.logical, sum) %>%
    tidyr::gather('celltype', 'specific', -c('family')) %>%
    dplyr::mutate(celltype=factor(celltype, levels=celltypes.rep)) %>%
    dplyr::select(family,specific,celltype)

stopifnot( all(loccount.celltype$family == loccount.celltype.specific$family) )
stopifnot( all(loccount.celltype$celltype == loccount.celltype.specific$celltype) )

loccount.celltype.specific <- loccount.celltype.specific %>%
    mutate(shared=loccount.celltype$count - specific)

p.all <- loccount.celltype.specific %>%
    dplyr::select(family,celltype,specific,shared) %>%    
    tidyr::gather('key', 'count', 3:4) %>%
    dplyr::mutate(key=factor(key, levels=c('specific', 'shared'))) %>%    
    ggplot(aes(celltype, count, label=count)) +
        geom_col(aes(fill=celltype, alpha=key), color='#333333') +
        geom_text(size = 3, position = position_stack(vjust = 0.5), color='#010101') +   
        scale_fill_manual(values=celltype_colors, guide=FALSE) +
        scale_alpha_discrete(name='', labels=c('Specific', 'Shared'), range=c(1.0, 0.60), guide=FALSE) +    
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~family, ncol=1, scales='free_y', strip.position='right') +
        xlab('') + ylab('# HERV elements expressed')

if(makepdf) pdf('analysis/figures/loccount_specific_ALL.celltypeXfamily.pdf', width=4, height=120)
p.all
if(makepdf) dev.off()

p <- loccount.celltype.specific %>%
    dplyr::select(family,celltype,specific,shared) %>%    
    tidyr::gather('key', 'count', 3:4) %>%
    dplyr::mutate(key=factor(key, levels=c('specific', 'shared'))) %>%    
    dplyr::filter(family %in% selfams) %>%
    dplyr::mutate(family=factor(family, levels=selfams)) %>%
    ggplot(aes(celltype, count, label=count)) +
        geom_col(aes(fill=celltype, alpha=key), color='#333333') +
        geom_text(size = 3, position = position_stack(vjust = 0.5), color='#010101') +   
        scale_fill_manual(values=celltype_colors, guide=FALSE) +
        scale_alpha_discrete(name='', labels=c('Specific', 'Shared'), range=c(1.0, 0.60), guide=FALSE) +    
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~family, ncol=1, scales='free_y', strip.position='right') +
        xlab('') + ylab('# HERV elements expressed')

if(makepdf) pdf('analysis/figures/loccount_specific.celltypeXfamily.pdf', paper='letter', width=4, height=11)
p
if(makepdf) dev.off()

################################################################################
# Total CPM, celltype X family
################################################################################
cpm_fam <- cpm.telescope %>% data.frame %>%
    tibble::rownames_to_column('locus') %>%
    dplyr::mutate(family=annot.herv$family) %>%
    group_by(family) %>%
    summarize_if(is.numeric, sum) %>%
    tidyr::gather('sample', 'cpm', -c('family')) %>%
    mutate(
        celltype=samples[sample,]$celltype,
        rnaextract=samples[sample,]$rnaextract
    ) %>%
    filter(rnaextract %in% extracts & celltype %in% celltypes.rep)

p.all <- cpm_fam %>%
    ggplot(aes(x=celltype, y=cpm)) + 
    geom_boxplot(aes(fill=celltype), outlier.shape=NA) + 
    geom_jitter(aes(color=celltype), height=0, width=0.2) +
    scale_fill_manual(values=celltype_colors, guide=FALSE) +
    scale_color_manual(values=celltype_colors, guide=FALSE) +    
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~family, ncol=1, scales='free_y', strip.position='right') +
    xlab('') + ylab('CPM')

if(makepdf) pdf('analysis/figures/cpm_ALL.celltypeXfamily.pdf', width=4, height=120)
p.all
if(makepdf) dev.off()

p <- cpm_fam %>%
    filter(family %in% selfams) %>%
    dplyr::mutate(family=factor(family, levels=selfams)) %>%    
    ggplot(aes(x=celltype, y=cpm)) + 
        geom_boxplot(aes(fill=celltype), color='#333333', outlier.shape=NA) + 
        geom_jitter(height=0, width=0.1, color='#333333', size=1, shape=18) +
        scale_fill_manual(values=celltype_colors, guide=FALSE) +
        scale_color_manual(values=celltype_colors, guide=FALSE) + 
        scale_y_continuous(position = "right") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~family, ncol=1, scales='free_y', strip.position='left') +
        xlab('') + ylab('CPM')

if(makepdf) pdf('analysis/figures/cpm.celltypeXfamily.pdf', paper='letter', width=4, height=11)
p
if(makepdf) dev.off()

