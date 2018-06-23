#! /usr/bin/env Rscript
library(tidyverse)

################################################################################
# Load sample data
################################################################################

#--- Cell type data (cancer vs. normal, tissue, etc.) is provided in table
#    named metadata/cell_lines.tsv
cell_metadata <- read.table('metadata/cell_lines.tsv', header=T, sep='\t') %>%
    mutate(celltype=factor(celltype, levels=celltype))

# Possible subcellular fraction for complete data
subcellular_fraction_levels <- 
    c('cell', 'cytosol', 'nucleus', 'nucleolus', 'nucleoplasm', 'chromatin')

#--- Load set1 metadata
samples.long <- read.table('metadata/set1.tsv',
           header=T, sep='\t', stringsAsFactors = F, comment.char="") %>%
    mutate(
        sample = Sample_Name,
        celltype=factor(source_name, levels=cell_metadata$celltype),
        layout=factor(LibraryLayout, 
                      levels=c('PAIRED','SINGLE')),
        rnaextract=factor(rnaextract, 
                          levels=c('longPolyA', 'total', 'longNonPolyA')),
        subcellular_fraction=factor(subcellular_fraction, 
                                    levels=subcellular_fraction_levels),
        sequencing_center=factor(sequencing_center,
                                 levels=c('CSHL','CALTECH'))
    )

row.names(samples.long) <- samples.long$sample

samples <- samples.long %>%
    select(sample, celltype, layout, rnaextract, subcellular_fraction, 
           sequencing_center)

samples <- samples %>% dplyr::inner_join(cell_metadata, by='celltype')

row.names(samples) <- samples$sample

# Clean up factors (for set1 only)
cell_metadata <- cell_metadata %>% 
    filter(celltype %in% samples$celltype) %>% 
    droplevels

samples <- samples %>% droplevels

rm(subcellular_fraction_levels)
