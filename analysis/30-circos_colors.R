#! /usr/bin/env Rscript
source('analysis/09-colors.R')
circos_dir <- 'analysis/circos'

################################################################################
# Setup colors
################################################################################

outfile <- file.path(circos_dir, 'mycolors.conf')

x <- col2rgb(heatmap_colors)
cat(
    sprintf('heatcolor-min = %s', 
            paste(col2rgb(RColorBrewer::brewer.pal(9, 'YlOrRd')[2])[,1], collapse=',')
    ),
    sprintf('heatcolor-max = %s', 
            paste(col2rgb(RColorBrewer::brewer.pal(9, 'YlOrRd')[9])[,1], collapse=',')
    ),
    sapply(1:length(heatmap_colors), function(col) 
        sprintf("heatcolor-%d-seq-%d = %d,%d,%d", 
                length(heatmap_colors), col, x[1,col],x[2,col],x[3,col])
    ),
    sep='\n',
    file=outfile, append=FALSE
)

cat(
    sapply(names(celltype_colors), function(cur_type) {
        tmp <- col2rgb(celltype_colors[cur_type])
        sprintf("cell-color-%s = %d,%d,%d", 
                tolower(gsub('\\+', '', cur_type)), tmp[1,1], tmp[2,1], tmp[3,1]
        )
    }),
    sep='\n',
    file=outfile, append=TRUE
)

rm(x)
