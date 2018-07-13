#! /usr/bin/env Rscript
library(tidyverse)

################################################################################
# Load annotation data
################################################################################

# Load mapping between HERV families and nicknames. This is used for coloring
# and visualization.
herv_fam <- read.table('analysis/herv_families.tsv', 
                       sep='\t', header=T, stringsAsFactors = F)

# Load TSV file with one row for each HERV in the annotation. Telescope reports
# are put into the same order as this TSV, and nicknames are added to the table.
annot.herv <- read.table('refs/HERV_rmsk.hg38.v2.tsv',
                         sep='\t', header=T, stringsAsFactors = F) %>%
    dplyr::mutate(length=end - start + 1) %>%
    tidyr::separate(locus, c("family"), sep='_', extra='drop', remove=F) %>%
    dplyr::left_join(herv_fam, by='family') %>%
    dplyr::mutate(
        family=factor(family, levels=herv_fam$family),
        group=factor(group, levels=unique(herv_fam$group)),
        letter=factor(letter, levels=unique(herv_fam$letter))
    ) %>%
    dplyr::select(locus, chrom, start, end, length, family, group, letter, category)
