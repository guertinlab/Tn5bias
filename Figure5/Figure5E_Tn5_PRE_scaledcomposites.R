library(lattice)
source('composite_functions.R')
################################################################################
################################################################################
################################################################################
################################################################################

#
train = c('AR', 'ASCL1', 'ATF3', 'CEBPB', 'CLOCK', 'CTCF', 'DLX2', 'DUX4', 'EGR1', 'ELF1', 'ESR1', 'FOXA1', 'FOXK2', 'GLIS1', 'HSF1', 'JUN', 'LEF1', 'MEIS2', 'MLX', 'MYC', 'NR2F2', 'SPIB', 'SRY', 'STAT1', 'TEAD1', 'TGIF1')
test = c('E2F1', 'EBF1', 'FERD3L', 'FOS', 'GATA2', 'HOXC12', 'IRF1', 'MAX', 'MEF2A', 'MGA', 'NFATC3', 'POU3F1', 'PPARG', 'REST', 'RUNX1', 'Six3', 'SP1', 'USF1')



load('Figure5E_motiflengths.Rdata')

load('Figure5E_Tn5_PRE_NNNCNN_unscaled_composites.Rdata')
#trim composites to 20bp around center
for (i in 1:length(combined_compositelist)) {
  combined_compositelist[[i]] = combined_compositelist[[i]][which(combined_compositelist[[i]]$x >= -20.5 & combined_compositelist[[i]]$x <= 20.5),]
}


for (i in c(10,13,35,44)) {
  plot.composites(combined_compositelist[[i]], legend = TRUE, 
                  pdf_name = paste('Figure5E_Tn5_PRE_NNNCNN_comparison_', paste(names(combined_compositelist[i])), '_composite', sep = ''),
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = TRUE, Motiflen = Motiflen[i]
  )}


