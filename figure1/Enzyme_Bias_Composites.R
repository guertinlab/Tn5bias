library(bigWig)
library(zoo)
library(lattice)
library(data.table)
library(matrixStats)
library(gridExtra)
source('../scripts/composite_functions.R')
################################################################################
################################################################################
################################################################################
################################################################################

setwd('../git')
################################################################################
################################################################################
load('../git/Unscaled_enzyme_composites.Rdata')


#Make vector of each enzyme bias motif's length, with the rowname the vector's name
mlen <- c(7, 5, 6, 4, 19)
mlen <- as.data.frame(mlen)
rownames(mlen) <- levels(as.factor(unscaled_composite_agg$factor))



enzymebias_plot <- plot.composites(unscaled_composite_agg, legend = FALSE, 
                ylabel = 'Cut Frequency',
                xlabel = 'Distance from Motif Center',
                figwidth = 4, figheight = 12,
                motifline = TRUE,
                Motiflen = mlen,
                striplabel = TRUE,
                pdf_name = 'enzymebias_composites',
                layoutgrid = c(1,5)
                
)


