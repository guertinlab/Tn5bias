library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
########################### For only plus unscaled composites
#Load in plus coordinates for each TF Motif
Motifs <- list.files('./')
Motifs = Motifs[grep('plus_400k_fimo', Motifs)]
Motiflist <- vector('list', length(Motifs))
for (i in 1:length(Motifs)) {
  Motiflist[[i]] <- FIMO.to.BED(paste('./',Motifs[i], sep = ''))
}
names(Motiflist) = substr(Motifs, 1, nchar(Motifs)-14)
Motifs = names(Motiflist)
#Determine signal at each plus motif
BWs <- c('../Figure1/DNase_Naked_plus_unscaled.bigWig')
DNase_plus_unscaled_compositelist <- vector('list', length(Motifs))
for (i in 1:length(Motiflist)) {
  DNase_plus_unscaled_compositelist[[i]] = BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                           upstream = 100, downstream = 100, factor = Motifs[i],
                                           group = 'Unscaled', ATAC = FALSE)
}
names(DNase_plus_unscaled_compositelist) <- Motifs

#plot uncorrected values
for (i in 1:length(DNase_plus_unscaled_compositelist)) {
  plot.composites(DNase_plus_unscaled_compositelist[[i]], legend = FALSE, 
                  pdf_name = paste('DNase_plus_unscaled_400k_', 
                  paste(names(DNase_plus_unscaled_compositelist[i])), '_composite', sep = ''),
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = Motiflen[i],
                  y_axis_range = c(0, max(DNase_plus_unscaled_compositelist[[i]]$est)),
                  x_axis_range = -100:100, nXticks = 4,
                  figwidth = 4, figheight=4, nYticks = 2, nYaxisdigits = 2
  )}
#Save values
save(DNase_plus_unscaled_compositelist, file = 'DNase_plus_unscaled_compositelist.Rdata')

########################### For only minus unscaled composites
#Load in minus coordinates for each TF Motif
Motifs <- list.files('./')
Motifs = Motifs[grep('minus_400k_fimo', Motifs)]
Motiflist <- vector('list', length(Motifs))
for (i in 1:length(Motifs)) {
  Motiflist[[i]] <- FIMO.to.BED(paste('./',Motifs[i], sep = ''))
}
names(Motiflist) = substr(Motifs, 1, nchar(Motifs)-14)
Motifs = names(Motiflist)
#Determine signal at each minus motif
BWs <- c('../Figure1/DNase_Naked_minus_unscaled.bigWig')
DNase_minus_unscaled_compositelist <- vector('list', length(Motifs))
for (i in 1:length(Motiflist)) {
  DNase_minus_unscaled_compositelist[[i]] = BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                            upstream = 100, downstream = 100, factor = Motifs[i],
                                            group = 'Unscaled', ATAC = FALSE)
}
names(DNase_minus_unscaled_compositelist) <- Motifs

#plot uncorrected values
for (i in 1:length(DNase_minus_unscaled_compositelist)) {
  plot.composites(DNase_minus_unscaled_compositelist[[i]], legend = FALSE, 
                  pdf_name = paste('DNase_minus_unscaled_400k_', 
                  paste(names(DNase_minus_unscaled_compositelist[i])), '_composite', sep = ''),
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = Motiflen[i],
                  y_axis_range = c(0, max(DNase_plus_unscaled_compositelist[[i]]$est)),
                  x_axis_range = -100:100, nXticks = 4,
                  figwidth = 4, figheight=4, nYticks = 2, nYaxisdigits = 2
  )}
#Save values
save(DNase_minus_unscaled_compositelist, file = 'DNase_minus_unscaled_compositelist.Rdata')