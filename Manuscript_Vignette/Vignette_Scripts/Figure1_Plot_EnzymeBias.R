source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
###################################################################################
library(bigWig)
library(zoo)
library(lattice)
###########################Benzonase unsep/minus/plus
Motifs <- c('Benzonase_sepcat_bias_transfac_00005_plus_400k_fimo.txt',
            'Benzonase_sepcat_bias_transfac_00005_minus_400k_fimo.txt')
Motiflist <- vector('list', length(Motifs))

for (i in 1:length(Motifs)) {
  Motiflist[[i]] = FIMO.to.BED(Motifs[i])
}
Motiflist[[1]] = rbind(Motiflist[[1]], Motiflist[[2]])

names(Motiflist) <- 'Benzonase'

BWs <- c('mm39_liver_Benzonase_plus.bigWig', 'mm39_liver_Benzonase_minus.bigWig')

merged_benz <- BED.query.bigWig(Motiflist[[1]], BWs[1], BWs[2], 
                                upstream = 50, downstream = 50, factor = 'Benzonase',
                                group = 'Unscaled', ATAC = FALSE)

unscaled_composite_agg <- rbind(merged_benz)

###########################Cyanase unsep/minus/plus
Motifs <- c('Cyanase_sepcat_bias_transfac_0001_plus_400k_fimo.txt',
            'Cyanase_sepcat_bias_transfac_0001_minus_400k_fimo.txt')
Motiflist <- vector('list', length(Motifs))

for (i in 1:length(Motifs)) {
  Motiflist[[i]] = FIMO.to.BED(Motifs[i])
}
Motiflist[[1]] = rbind(Motiflist[[1]], Motiflist[[2]])

names(Motiflist) <- 'Cyanase'

BWs <- c('mm39_liver_Cyanase_plus.bigWig', 'mm39_liver_Cyanase_minus.bigWig')

merged_cyan <- BED.query.bigWig(Motiflist[[1]], BWs[1], BWs[2], 
                                upstream = 50, downstream = 50, factor = 'Cyanase',
                                group = 'Unscaled', ATAC = FALSE)


unscaled_composite_agg <- rbind(unscaled_composite_agg, merged_cyan)

###########################DNase unsep/minus/plus
Motifs <- c('DNase_sepcat_bias_transfac_0001_plus_400k_fimo.txt',
            'DNase_sepcat_bias_transfac_0001_minus_400k_fimo.txt')
Motiflist <- vector('list', length(Motifs))

for (i in 1:length(Motifs)) {
  Motiflist[[i]] = FIMO.to.BED(Motifs[i])
}
Motiflist[[1]] = rbind(Motiflist[[1]], Motiflist[[2]])
names(Motiflist) <- 'DNase'

BWs <- c('DNase_Naked_plus_unscaled.bigWig', 'DNase_Naked_minus_unscaled.bigWig')

merged_DNase <- BED.query.bigWig(Motiflist[[1]], BWs[1], BWs[2], 
                                 upstream = 50, downstream = 50, factor = 'DNase',
                                 group = 'Unscaled', ATAC = FALSE)


unscaled_composite_agg <- rbind(unscaled_composite_agg, merged_DNase)

###########################MNase unsep/minus/plus
Motifs <- c('MNase_sepcat_bias_transfac_0001_plus_400k_fimo.txt',
            'MNase_sepcat_bias_transfac_0001_minus_400k_fimo.txt')
Motiflist <- vector('list', length(Motifs))

for (i in 1:length(Motifs)) {
  Motiflist[[i]] = FIMO.to.BED(Motifs[i])
}
Motiflist[[1]] = rbind(Motiflist[[1]], Motiflist[[2]])
names(Motiflist) <- 'MNase'

BWs <- c('mm39_liver_MNase_plus.bigWig', 'mm39_liver_MNase_minus.bigWig')

merged_MNase <- BED.query.bigWig(Motiflist[[1]], BWs[1], BWs[2], 
                                 upstream = 50, downstream = 50, factor = 'MNase',
                                 group = 'Unscaled', ATAC = FALSE)


unscaled_composite_agg <- rbind(unscaled_composite_agg, merged_MNase)

###########################Tn5 unsep/minus/plus
Motifs <- c('Tn5_sepcat_bias_transfac_00005_plus_400k_fimo.txt',
            'Tn5_sepcat_bias_transfac_00005_minus_400k_fimo.txt')
Motiflist <- vector('list', length(Motifs))

for (i in 1:length(Motifs)) {
  Motiflist[[i]] = FIMO.to.BED(Motifs[i])
}
Motiflist[[1]] = rbind(Motiflist[[1]], Motiflist[[2]])
names(Motiflist) <- 'Tn5'

BWs <- c('C1_gDNA_rep1_plus.bigWig', 'C1_gDNA_rep1_minus.bigWig')

merged_Tn5 <- BED.query.bigWig(Motiflist[[1]], BWs[1], BWs[2], 
                               upstream = 50, downstream = 50, factor = 'Tn5',
                               group = 'Unscaled', ATAC = TRUE)

unscaled_composite_agg <- rbind(unscaled_composite_agg, merged_Tn5)

#Determine each enzyme's read depth:
benzonase_read_depth = fread('mm39_liver_Benzonase_not_scaled.bed')
benzonase_read_depth = sum(benzonase_read_depth$V5)
cyanase_read_depth = fread('mm39_liver_Cyanase_not_scaled.bed')
cyanase_read_depth = sum(cyanase_read_depth$V5)
MNase_read_depth = fread('mm39_liver_MNase_not_scaled.bed')
MNase_read_depth = sum(MNase_read_depth$V5)
DNase_read_depth = fread('DNase_Naked_unscaled_not_scaled.bed')
DNase_read_depth = sum(DNase_read_depth$V5)
Tn5_read_depth = fread('C1_gDNA_rep1_not_scaled.bed')
Tn5_read_depth = sum(Tn5_read_depth$V5)


#Scale all values by their read depth
unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'Benzonase')] = 
  unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'Benzonase')]/benzonase_read_depth

unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'Cyanase')] = 
  unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'Cyanase')]/cyanase_read_depth

unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'MNase')] = 
  unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'MNase')]/MNase_read_depth

unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'DNase')] = 
  unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'DNase')]/DNase_read_depth

unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'Tn5')] = 
  unscaled_composite_agg$est[which(unscaled_composite_agg$factor == 'Tn5')]/Tn5_read_depth

#Normalize all plots to greatest value
unscaled_composite_agg$est = unscaled_composite_agg$est/max(unscaled_composite_agg$est)

###Plot these all together
enzymebiasplot = plot.composites(unscaled_composite_agg, x_axis_range=-15:15,
                          legend = FALSE, 
                          ylabel = 'Cut Frequency',
                          xlabel = 'Distance from Motif Center',
                          figwidth = 4, figheight = 12,
                          motifline = FALSE,
                          Motiflen = mlen,
                          striplabel = TRUE,
                          pdf_name = 'unscaled_enzymebias_composites',
                          indexlist = list(c(5,3,2,1,4)),
                          layoutgrid = c(1,5),
                          y_axis_range = 0:1, y_axis = TRUE, Y_axis_ticks = seq(0,1,0.5),
)

