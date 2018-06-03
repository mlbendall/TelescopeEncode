library(tidyverse)

cell_metadata <- read.table('metadata/cell_lines.tsv', header=T, sep='\t') %>%
    mutate(celltype=factor(celltype, levels=celltype))


subcellular_fraction_levels <- 
    c('cell', 'cytosol', 'nucleus', 'nucleolus', 'nucleoplasm', 'chromatin')


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

rm(subcellular_fraction_levels)
