#! /usr/bin/env Rscript
source('analysis/07-expressed_cpm.R')

library(tidyverse)
library(cowplot)

makepdf <- TRUE

################################################################################
# Make HERV family profiles
################################################################################
prop_cutoff <- 0.05

# Mean CPM (locus x celltype)
cpm_celltype <- lapply(celltypes.rep, function(ctype){
    tmp <- cpm.telescope[ , samples$rnaextract %in% extracts & samples$celltype == ctype, drop=FALSE]
    data.frame(rowMeans(tmp))
}) %>% bind_cols()
colnames(cpm_celltype) <- celltypes.rep
row.names(cpm_celltype) <- row.names(cpm.telescope)

# Proportion of retrotranscriptome (locus x celltype)
prop_celltype <- data.frame(t(t(cpm_celltype) / colSums(cpm_celltype)))

# Proportion by family (family x celltype)
prop_celltype_fam <- prop_celltype %>%
    tibble::rownames_to_column('locus') %>%
    dplyr::mutate(
        hervfam = annot.herv$family,
    )  %>%
    group_by(hervfam) %>%
    summarize_if(is.numeric, sum) %>%
    data.frame

row.names(prop_celltype_fam) <- prop_celltype_fam$hervfam
prop_celltype_fam$hervfam <- NULL

bigfams <- prop_celltype_fam[rowSums(prop_celltype_fam > prop_cutoff ) > 0,]
rem <- prop_celltype_fam[rowSums(prop_celltype_fam > prop_cutoff) == 0,]

# Assign groups to families
rem$hervgrp <- sapply(row.names(rem), function(x) herv_fam[herv_fam$family==x,]$group)

prop_celltype_grp <- rem %>%
    group_by(hervgrp) %>%
    summarize_if(is.numeric, sum) %>%
    data.frame

row.names(prop_celltype_grp) <- prop_celltype_grp$hervgrp
prop_celltype_grp$hervgrp <- NULL

biggrps <- prop_celltype_grp[rowSums(prop_celltype_grp > prop_cutoff ) > 0,]
rem2 <- prop_celltype_grp[rowSums(prop_celltype_grp > prop_cutoff ) == 0,]

stopifnot(all((colSums(bigfams) + colSums(biggrps) + colSums(rem2)) == 1))

row.names(biggrps) <- paste0(row.names(biggrps), '.grp')

famlevels <- c(
    row.names(bigfams)[order(matrixStats::rowMaxs(as.matrix(bigfams)), decreasing = T)],
    row.names(biggrps)[order(matrixStats::rowMaxs(as.matrix(biggrps)), decreasing = T)],
    "other"
)
hervfam_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(length(famlevels))
names(hervfam_colors) <- famlevels

prop_celltype_combined <- rbind(bigfams, biggrps, other=colSums(rem2)) %>% data.frame
colnames(prop_celltype_combined) <- gsub('\\.', '-', colnames(prop_celltype_combined))

p <- prop_celltype_combined %>%
    tibble::rownames_to_column('family') %>%
    tidyr::gather('celltype', 'prop', -c('family')) %>%
    dplyr::mutate(
        family=factor(family, levels=famlevels),
        celltype=factor(celltype, levels=celltypes.rep)
    ) %>%
    ggplot(aes(x=celltype)) +
    geom_col(aes(y=prop, fill=family), 
             position=position_stack(reverse=TRUE)
    ) +
    scale_fill_manual(values=hervfam_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('') + ylab('Proportion of HERV expression')

if(makepdf) pdf('analysis/figures/herv_family_profiles.pdf', paper='USr', width=10, height=5)
p
if(makepdf) dev.off()

################################################################################
# Make top HERV for each celltype
################################################################################
nlocs <- 10

plots.list <- lapply(celltypes.rep, function(ctype){
    meanexp <- data.frame(
        locus=annot.herv$locus,
        family=annot.herv$family,
        cpm=cpm_celltype[,ctype]
    )
    tmp <- cpm.telescope[ ,samples$rnaextract %in% extracts & samples$celltype == ctype, drop=FALSE]
    meanexp$sd <- sapply(1:nrow(tmp), function(x) sd(tmp[x,])) 
    meanexp$se <- meanexp$sd / sqrt(ncol(tmp))  
    
    if(cell_metadata[cell_metadata$celltype==ctype,]$sex == 'F')
        meanexp[annot.herv$chrom == 'chrY', c(3,4,5)] <- 0
    toplocs <- order(meanexp$cpm, decreasing = TRUE)[1:nlocs]

    p1 <- meanexp[toplocs,] %>%
        dplyr::mutate(locus=factor(locus, levels=locus)) %>%
        ggplot(aes(x=locus, y=cpm)) +
            geom_col(color='#333333', fill='#CCCCCC') +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
            xlab('') + ylab('CPM') + ggtitle(paste0(ctype))
    
    p2 <- p1 + geom_errorbar(aes(ymin=cpm, ymax=cpm+se), width=0.2, color='#333333')    
    
    allpoints <- tmp[toplocs,] %>%
        data.frame %>%
        tibble::rownames_to_column('locus') %>%
        tidyr::gather('sample', 'cpm', -c(locus)) %>%
        mutate(locus=factor(locus, levels=annot.herv[toplocs,]$locus))
    
    p3 <- p2 + geom_point(data=allpoints, aes(x=locus, y=cpm))
    
    p2
})
names(plots.list) <- celltypes.rep

if(makepdf) pdf('analysis/figures/topherv_each.pdf', paper='USr', width=4, height=2.5)
for(p in plots.list){
    plot(p)
}
if(makepdf) dev.off()
    

# Create plotgrid with titles and axis removed
plotgrid <- plot_grid(plotlist=lapply(plots.list, function(p) {
    p + theme(
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8))
    }),
    ncol=4,nrow=2,align='hv'
)

save_plot(paste0('analysis/figures/topherv_combined.pdf'),
          plotgrid, nrow=2, base_height=2.5, base_width=16)

