library(ggplot2)
source('composite_functions.R')
################################################################################
################################################################################
################################################################################
################################################################################
train = c('AR', 'ASCL1', 'ATF3', 'CEBPB', 'CLOCK', 'CTCF', 'DLX2', 'DUX4', 'EGR1', 'ELF1', 'ESR1', 'FOXA1', 'FOXK2', 'GLIS1', 'HSF1', 'JUN', 'LEF1', 'MEIS2', 'MLX', 'MYC', 'NR2F2', 'SPIB', 'SRY', 'STAT1', 'TEAD1', 'TGIF1')
test = c('E2F1', 'EBF1', 'FERD3L', 'FOS', 'GATA2', 'HOXC12', 'IRF1', 'MAX', 'MEF2A', 'MGA', 'NFATC3', 'POU3F1', 'PPARG', 'REST', 'RUNX1', 'Six3', 'SP1', 'USF1')



load('Figure3C_DNase_PRE_NNNCNN_unscaled_compositelist.Rdata')
load('Figure3C_PRE_TF_Motiflen.Rdata')
#trim composites to 20bp around center
for (i in 1:length(combined_compositelist)) {
  combined_compositelist[[i]] = combined_compositelist[[i]][which(combined_compositelist[[i]]$x >= -20.5 & combined_compositelist[[i]]$x <= 20.5),]
}

#Separate out the TFs to be plotted and their motif lengths
figure3C_plot = do.call(rbind, combined_compositelist[c(13, 35, 44)])
mlen <- Motiflen[c(13, 35, 44)]
mlen <- as.data.frame(mlen)
rownames(mlen) <- levels(as.factor(figure3C_plot$factor))



plot.composites(figure3C_plot, legend = TRUE, 
                  pdf_name = 'Figure3C_DNasePN_PRE_NNNCNN_comparison',
                  figwidth = 8, figheight = 4,
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = TRUE, Motiflen = mlen, layoutgrid = c(3,1))

