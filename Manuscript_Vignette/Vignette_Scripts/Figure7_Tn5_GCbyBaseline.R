source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
library(ggplot2)

################################################################################
load('Tn5_unscaled_compositelist.Rdata')
load('../Figure5/Motiflen.Rdata')
#Remove single nucleotide values
unscaled_composites = do.call(rbind, Tn5_unscaled_compositelist)

baseline_frame = NULL
baseline_store = NULL
for (i in 1:length(Motiflen)) {
  baseline_store = unscaled_composites[which(unscaled_composites$factor == names(Motiflen)[i] &
                   unscaled_composites$x > ((Motiflen[i]/2)+10) | unscaled_composites$factor ==
                   names(Motiflen)[i] & unscaled_composites$x < -((Motiflen[i]/2)+10)),]
  baseline_frame = rbind(baseline_frame, baseline_store)
}
unscaled_composites = baseline_frame

unscaled_composites$factor = as.factor(unscaled_composites$factor)

#Determine baseline mean for each TF
baseline_mean = aggregate(unscaled_composites[,2], list(unscaled_composites$factor), mean)

#Plot unscaled and PRE output GC content
Motifs = list.files('../Figure4/')
Motifs = Motifs[grep('meme',Motifs)]
Motifs = Motifs[-c(8)]

Motifs = substr(Motifs, 1, nchar(Motifs)-5)
GCvec <- NULL
for(i in 1:length(Motifs)) {
  memefile <- paste0("../Figure4/",Motifs[i], ".meme")
  minmeme <- read.csv(memefile, skip = 11, header = FALSE, sep = '')
  minmeme <- minmeme[1:nrow(minmeme)-1,]
  minmeme[,2] <- as.numeric(minmeme[,2])
  minmeme[,1] <- as.numeric(minmeme[,1])
  GCcon <- sum(minmeme[,c(2,3)])/sum(minmeme)
  GCvec[i] <- GCcon
}
baseline_mean[,3] = GCvec
colnames(baseline_mean) = c('TF', 'BaselineAvg', 'GC')
pdf("Figure7B_Tn5_unscaled_BaselineVsGCcon.pdf", width=6, height=6)
ggplot(baseline_mean, aes(x = BaselineAvg, y = GC, label = TF)) +
  geom_point(stat = 'identity') +
  xlim(0, 0.016) +
  ylim(0, 0.8) +
  theme_classic() + xlab('Baseline Signal') + ylab('Motif GC%') + 

  theme(axis.text=element_text(size=20, color = 'black'),
        axis.title=element_text(size=30, color = 'black'),
        axis.text.y=element_text(hjust = 0.4),
        axis.text.x=element_text(hjust = 0.60)) 
dev.off()

