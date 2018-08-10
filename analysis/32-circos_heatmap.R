#! /usr/bin/env Rscript
source('analysis/07-expressed_cpm.R')
source('analysis/09-colors.R')

library(tidyverse)

circos_dir <- 'analysis/circos'
ok_chroms <- c(sprintf('chr%d', 1:22), 'chrX', 'chrY')

################################################################################
# Setup celltype radius
################################################################################
r1_start <- 0.90
ctype_size <- 0.06
ctype_space <- 0.005

ctype_r1 <- seq(r1_start, length.out = length(celltypes.rep), by=-ctype_size)
ctype_r0 <- ctype_r1 - (ctype_size - ctype_space)
names(ctype_r1) <- celltypes.rep
names(ctype_r0) <- celltypes.rep

#--- Generate rings.conf
outfile <- file.path(circos_dir, 'rings.conf')
cat('<plot>',
    sprintf('r1=%0.03fr', r1_start + 0.004),
    sprintf('r0=%0.03fr', r1_start + 0.002),
    '<backgrounds>', '<background>',
    'color = greys-3-seq-3',
    '</background>','</backgrounds>',
    '</plot>',
    sep='\n',
    file=outfile, append=FALSE
)

for (cur_type in celltypes.rep) {
    cat('<plot>',
        sprintf('r1=%0.03fr', ctype_r0[cur_type] - 0.002),
        sprintf('r0=%0.03fr', ctype_r0[cur_type] - 0.004),
        '<backgrounds>', '<background>',
        # sprintf('color = conf(colors, cell-color-%s),0.75', tolower(gsub('\\+','', cur_type))),        
        'color = vlgrey',
        '</background>','</backgrounds>',
        '</plot>',
        sep='\n',
        file=outfile, append=TRUE
    )
}

outfile <- file.path(circos_dir, 'celltype_backgrounds.conf')
cat('#\n', file=outfile, append=FALSE)
for (cur_type in celltypes.rep) {
    cat('<plot>',
        sprintf('r1=%0.03fr', ctype_r1[cur_type] + 0.003),
        sprintf('r0=%0.03fr', ctype_r0[cur_type] - 0.003),
        'z = -10',
        '<backgrounds>', '<background>',
        sprintf('color = conf(colors, cell-color-%s),0.95', tolower(gsub('\\+','', cur_type))),
        '</background>','</backgrounds>',
        '</plot>',
        sep='\n',
        file=outfile, append=TRUE
    )
}

################################################################################
# HERV expression heatmap - by sample
################################################################################
sfilt <- (samples$rnaextract %in% extracts) & (samples$celltype %in% celltypes.rep)

lapply(samples[sfilt,]$sample, function(s) {
    cur_type <- samples[s,'celltype']
    rel_f <- file.path('data', paste(s, 'heatmap', 'txt', sep='.'))
    outfile <- file.path(circos_dir, rel_f)

    tmp <- annot.herv[,c('chrom','start','end','locus')]
    tmp$value <- log2(1 + cpm.telescope[,s])
    tmp <- tmp[expressed.celltype[ ,cur_type], ]
    tmp <- tmp[tmp$chrom %in% ok_chroms,]
    if(samples[s, 'sex'] == 'F') tmp <- tmp[tmp$chrom %in% ok_chroms[1:23], ]

    tmp <- tmp[tmp$value > 0, ]
    tmp %>%
        mutate(
            chrom=gsub('^chr', 'hs', chrom),
            opts=sprintf('locus=%s', locus)
        ) %>%
        select(chrom, start, end, value, opts) %>%
        write.table(file=outfile, sep='\t', row.names=F, col.names=F, quote=F)
    data.frame(sample=s, 
               celltype=cur_type, 
               filename=rel_f, 
               stringsAsFactors = F
               )
}) %>% bind_rows() -> filelist

allmin <- min(log2(1 + cpm.telescope[expressed.celltype.any , sfilt]))
allmax <- max(log2(1 + cpm.telescope[expressed.celltype.any , sfilt]))

sprintf("Minimum expression value: %0.5f", allmin)
sprintf("Maximum expression value: %0.5f", allmax)

filelist <- filelist %>% 
    dplyr::mutate(celltype=factor(celltype, levels=celltypes.rep)) %>%
    dplyr::arrange(celltype) %>% 
    dplyr::mutate(tmp=1) %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(repnum=cumsum(tmp), outof=sum(tmp)) %>%
    dplyr::select(-tmp) %>%
    mutate(r1 = ctype_r1[celltype] - ((repnum - 1) * ((ctype_size - ctype_space) / outof))) %>%
    mutate(r0 = r1 - ((ctype_size - ctype_space) / outof))

#--- Generate hm_sample.conf
cat(
    'type = heatmap',
    sprintf('color = heatcolor-min,%s,heatcolor-max',
            paste(sprintf('heatcolor-%d-seq-%d', length(heatmap_colors), 1:length(heatmap_colors)), collapse=',')
    ),
    'min = 0.0000',
    'max = 10.0000',
    'color_mapping = 2',
    'scale_log_base = 1',
    '',
    sep='\n',
    file=file.path(circos_dir, 'hm_sample.conf')
)

for(i in 1:nrow(filelist)) {
    row <- filelist[i,]
    cat(sprintf('<plot>\nfile=%s\nr1=%0.03fr\nr0=%0.03fr\n</plot>\n', 
                row$filename, row$r1, row$r0 ),
        sep='\n',
        file=file.path(circos_dir, 'hm_sample.conf'),
        append=TRUE
    )
}
