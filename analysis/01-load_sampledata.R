library(tidyverse)

celltype_levels <- c('H1-hESC', 
                     'GM12878',
                     'K562',
                     'HeLa-S3',
                     'HepG2',
                     'HUVEC',
                     'SK-N-SH',
                     'IMR90',
                     'A549',
                     'MCF-7',
                     'CD20+',
                     'Monocytes-CD14+',
                     'NHEK'
)
layout_levels <- c('PAIRED',
                   'SINGLE'
)
rnaextract_levels <- c('longPolyA', 
                       'total', 
                       'longNonPolyA'
)
localization_levels <- c('cell', 
                         'cytosol', 
                         'nucleus',
                         'nucleolus',                         
                         'nucleoplasm',
                         'chromatin'
)


samples <- read.table('metadata/set1.txt',
           header=T, sep='\t', stringsAsFactors = F, comment.char="") %>%
    mutate(
        sample_id = Sample_Name,
        celltype=factor(source_name, levels=celltype_levels),
        layout=factor(LibraryLayout, levels=layout_levels),
        rnaextract=factor(rnaextract, levels=rnaextract_levels),
        localization=factor(localization, levels=localization_levels)
    ) %>%
    select(sample_id, celltype, layout, rnaextract, localization)

row.names(samples) <- samples$sample_id
