#! /usr/bin/env Rscript

plotpal <- FALSE

# Function for plotting colors side-by-side
pal <- function(col, border = "light gray", ...){
    n <- length(col)
    plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "", ...)
    rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

################################################################################
# Palettes
################################################################################
colors.set1 <- RColorBrewer::brewer.pal(9, 'Set1')
names(colors.set1) <- c('red', 'blue', 'green', 'purple', 'orange',
                        'yellow', 'brown', 'pink', 'gray')

colors.pastel1 <- RColorBrewer::brewer.pal(9, 'Pastel1')
names(colors.pastel1) <- c('red', 'blue', 'green', 'purple', 'orange',
                           'yellow', 'brown', 'pink', 'gray')

cnames <- c('teal', 'vermillion', 'purple', 'magenta', 'green', 
            'gold', 'brown', 'gray')
colors.2 <- c(
    RColorBrewer::brewer.pal(8, 'Pastel2'),
    RColorBrewer::brewer.pal(8, 'Set2'),
    RColorBrewer::brewer.pal(8, 'Dark2')
)
names(colors.2) <- c(paste0('l', cnames), cnames, paste0('d', cnames))

colors.paired <- RColorBrewer::brewer.pal(12, 'Paired')
names(colors.paired) <- c('lblue', 'blue', 'lgreen', 'green', 'lred',
                          'red', 'lorange', 'orange', 'lpurple', 'purple',
                          'lbrown', 'brown')

colors.accent <- RColorBrewer::brewer.pal(8, 'Accent')
names(colors.accent) <- c('mint', 'lavender', 'peach', 'lemon', 'blue',
                          'magenta', 'brown', 'gray')

################################################################################
# Celltypes
################################################################################
celltype_colors <- c(
    "H1-hESC"          = colors.2[['dmagenta']],

    "HeLa-S3"          = colors.2[['teal']],
    "HepG2"            = colors.2[['purple']],
    "SK-N-SH"          = colors.accent[['mint']],
    "A549"             = colors.2[['green']],
    "MCF-7"            = colors.accent[['blue']],

    "HUVEC"            = colors.2[['dgreen']],
    "IMR90"            = colors.2[['dteal']],
    "NHEK"             = colors.2[['dpurple']],
    
    "GM12878"          = colors.2[['dvermillion']],
    "K562"             = colors.2[['vermillion']],
    "CD20+"            = colors.set1[['red']],
    "CD14+"  = colors.set1[['orange']]
)
if(plotpal) pal(celltype_colors)

karyotype_colors <- c(
    "NORMAL" = "#FFFFFF00",
    "CANCER" = "#333333"
)

layout_colors <- c(
    'SINGLE' = "#FFFFFF00",
    'PAIRED' = "#333333"
)

center_colors <- c(
    'CSHL' = "#333333",
    'CALTECH' = "#FFFFFF00"
)

################################################################################
# Heatmaps
################################################################################
# Bright colors to anchor heatmap
# Yellow (Gold) -> Orange -> Red
heatmap_anchors <- c(
    RColorBrewer::brewer.pal(8, 'Set2')[6],
    RColorBrewer::brewer.pal(9, 'Set1')[5],
    RColorBrewer::brewer.pal(9, 'YlOrRd')[8]
)
if(plotpal) pal(heatmap_anchors)

heatmap_colors <- colorRampPalette(heatmap_anchors)(50)
if(plotpal) pal(heatmap_colors)

rm(plotpal, pal, heatmap_anchors, colors.set1, colors.pastel1, colors.2, 
   colors.paired, colors.accent)

