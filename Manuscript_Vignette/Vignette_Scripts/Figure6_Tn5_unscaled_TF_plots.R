library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
########################### For only plus unscaled composites
#Load in plus coordinates for each TF Motif
Motifs <- list.files('../Figure5')
Motifs = Motifs[grep('plus_400k_fimo', Motifs)]
Motiflist <- vector('list', length(Motifs))
for (i in 1:length(Motifs)) {
  Motiflist[[i]] <- FIMO.to.BED(paste('../Figure5/',Motifs[i], sep = ''))
}
names(Motiflist) = substr(Motifs, 1, nchar(Motifs)-14)
Motifs = names(Motiflist)
#Determine signal at each plus motif
BWs <- c('../Figure1/C1_gDNA_rep1_plus.bigWig')
Tn5_plus_unscaled_compositelist <- vector('list', length(Motifs))
for (i in 1:length(Motiflist)) {
  Tn5_plus_unscaled_compositelist[[i]] = BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                           upstream = 100, downstream = 100, factor = Motifs[i],
                                           group = 'Unscaled', ATAC = TRUE)
}
names(Tn5_plus_unscaled_compositelist) <- Motifs

#plot uncorrected values
for (i in 1:length(Tn5_plus_unscaled_compositelist)) {
  plot.composites(Tn5_plus_unscaled_compositelist[[i]], legend = FALSE, 
                  pdf_name = paste('Tn5_plus_unscaled_400k_', 
                  paste(names(Tn5_plus_unscaled_compositelist[i])), '_composite', sep = ''),
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = Motiflen[i],
                  y_axis_range = c(0, max(Tn5_plus_unscaled_compositelist[[i]]$est)),
                  x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50),
                  figwidth = 4, figheight=4, nYticks = 2, nYaxisdigits = 2
  )}
#Save values
save(Tn5_plus_unscaled_compositelist, file = 'Tn5_plus_unscaled_compositelist.Rdata')

########################### For only minus unscaled composites
#Load in minus coordinates for each TF Motif
Motifs <- list.files('../Figure5')
Motifs = Motifs[grep('minus_400k_fimo', Motifs)]
Motiflist <- vector('list', length(Motifs))
for (i in 1:length(Motifs)) {
  Motiflist[[i]] <- FIMO.to.BED(paste('../Figure5/',Motifs[i], sep = ''))
}
names(Motiflist) = substr(Motifs, 1, nchar(Motifs)-14)
Motifs = names(Motiflist)
#Determine signal at each minus motif
BWs <- c('../Figure1/C1_gDNA_rep1_minus.bigWig')
Tn5_minus_unscaled_compositelist <- vector('list', length(Motifs))
for (i in 1:length(Motiflist)) {
  Tn5_minus_unscaled_compositelist[[i]] = BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                            upstream = 100, downstream = 100, factor = Motifs[i],
                                            group = 'Unscaled', ATAC = TRUE)
}
names(Tn5_minus_unscaled_compositelist) <- Motifs

#plot uncorrected values
for (i in 1:length(Tn5_minus_unscaled_compositelist)) {
  plot.composites(Tn5_minus_unscaled_compositelist[[i]], legend = FALSE, 
                  pdf_name = paste('Tn5_minus_unscaled_400k_', 
                  paste(names(Tn5_minus_unscaled_compositelist[i])), '_composite', sep = ''),
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = Motiflen[i],
                  y_axis_range = c(0, max(Tn5_plus_unscaled_compositelist[[i]]$est)),
                  x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50),
                  figwidth = 4, figheight=4, nYticks = 2, nYaxisdigits = 2
  )}
#Save values
save(Tn5_minus_unscaled_compositelist, file = 'Tn5_minus_unscaled_compositelist.Rdata')
