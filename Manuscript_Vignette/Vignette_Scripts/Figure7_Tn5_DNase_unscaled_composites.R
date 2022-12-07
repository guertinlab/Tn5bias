library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
library(bigWig)
library(lattice)
library(grid)
###########################
#Load in genomic locations determined using FIMO
Motifs <- list.files('../Figure5')
Motifs = Motifs[grep('_fimo', Motifs)]
Motiflist <- vector('list', length(Motifs))
for (i in 1:length(Motifs)) {
  Motiflist[[i]] <- FIMO.to.BED(paste('../Figure5/',Motifs[i], sep = ''))
}
names(Motiflist) = Motifs
#Combine the minus and plus FIMO motifs
minus_Motifs = Motiflist[seq(1, 88, 2)]
plus_Motifs = Motiflist[-c(seq(1, 88, 2))]
Motiflist = vector('list', length(plus_Motifs))
for (i in 1:length(plus_Motifs)) {
  Motiflist[[i]] = rbind(plus_Motifs[[i]], minus_Motifs[[i]])
}
names(Motiflist) = substr(names(plus_Motifs), 1 , nchar(names(plus_Motifs))-22)
rm(plus_Motifs, minus_Motifs)

#Make composite plots of unscaled values 
BWs <- c('../Figure1/C1_gDNA_rep1.bigWig')
compositelist <- vector('list', length(Motiflist))
for (i in 1:length(Motiflist)) {
  compositelist[[i]] = BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                        upstream = 100, downstream = 100,
                                        factor = names(Motiflist[i]),
                                        group = 'Unscaled', ATAC = TRUE)
}
names(compositelist) = names(Motiflist)
Tn5_unscaled_compositelist = compositelist
save(Tn5_unscaled_compositelist, file = 'Tn5_unscaled_compositelist.Rdata')
