library(data.table)
library(ggplot2)
library(seqLogo)
library(scclust)
library(RColorBrewer)
library(ggrepel)
library(ggnewscale)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
################################################################################
train = c( 'ASCL1', 'ATF3', 'CLOCK', 'CTCF',
           'DLX2', 'EBF1', 'EGR1', 'ELF1',
           'FOS', 'FOXA1', 'FOXK2', 'GATA2',
           'GLIS1', 'IRF1', 'JUN', 'LEF1',
           'MEIS2', 'MLX', 'MYC', 'PPARG', 
           'RUNX1', 'SPIB', 'SRY', 'STAT1', 'TGIF1')
test = c('AR', 'CEBPB', 'E2F1', 'ESR1',
         'FERD3L', 'HOXC12', 'HSF1', 'MAX',
         'MEF2A', 'MGA', 'NFATC3', 'NR2F2',
         'POU3F1', 'REST', 'Six3', 'SP1',
         'USF1', 'TEAD1')


#Determine GC and IC of each transcription factor
Memes <- list.files('./')
Memes = Memes[grep('.meme', Memes)]

GCvec = NULL
ICvec = NULL
for(i in 1:length(Memes)) {
  memefile <- paste0("./",Memes[i])
  minmeme <- read.csv(memefile, skip = 11, header = FALSE, sep = '')
  minmeme <- minmeme[1:nrow(minmeme)-1,]
  minmeme[,2] <- as.numeric(minmeme[,2])
  minmeme[,1] <- as.numeric(minmeme[,1])
  GCcon <- sum(minmeme[,c(2,3)])/sum(minmeme)
  GCvec[i] <- GCcon
  minmeme = t(minmeme)
  rownames(minmeme) = c('A', 'C', 'G', 'T')
  pwmMeme = makePWM(minmeme)
  ICcon = sum(ic(pwmMeme))
  ICvec[i] = ICcon
}

#Combine IC and GC vectors and add TF name column
GCbyIC = as.data.frame(cbind(GCvec, ICvec))
GCbyIC[,3] = substr(Memes, 1, nchar(Memes)- 5)
#Add column names
colnames(GCbyIC) = c('GC', 'IC', 'TF')
#z-standardize GC values
GCmean = mean(GCbyIC$GC)
GCstddv = sd(GCbyIC$GC)
GCbyICstd = as.data.frame((GCbyIC$GC-GCmean)/GCstddv)
#z-standardize IC values
ICmean = mean(GCbyIC$IC)
ICstddv = sd(GCbyIC$IC)
GCbyICstd[,2] = (GCbyIC$IC-ICmean)/ICstddv
GCbyICstd[,3] = GCbyIC$TF
#Add column names
colnames(GCbyICstd) = c('zstdGC', 'zstdIC', 'TF')

#Compute euclidean distances between points
GCbyIC_dist = distances(GCbyICstd,
                        id_variable = "TF",
                        dist_variables = c("zstdGC", "zstdIC"))
#Implement size-constrained clustering
set.seed(42)
GCbyIC_scclust = sc_clustering(GCbyIC_dist, 2)
#Refine clustering using hierarchical clustering
GCbyIC_clust = hierarchical_clustering(GCbyIC_dist, 2, batch_assign = TRUE,
                                       existing_clustering = GCbyIC_scclust)
#Add clusters to dataframe
GCbyICstd[,4] = GCbyIC_clust
colnames(GCbyICstd)[4] = 'Cluster'
GCbyICstd[,4] = as.character(GCbyICstd[,4])
#Remove DUX4
GCbyICstd = GCbyICstd[-c(8),]
#Assign test and training labels to each TF
GCbyICstd$group = ''
for (i in 1:nrow(GCbyICstd)) {
  if (GCbyICstd[i,3] %in% train) {
    GCbyICstd[i,5] = 'Train'
  } else {
    GCbyICstd[i,5] = 'Test'
  } 
}
#Assign colors to test and training sets and make a color pallette for clusters
GCbyICstd[,6] = 'red'
GCbyICstd[which(GCbyICstd[,5] == 'Train'),6] = 'blue'
colnames(GCbyICstd)[6] = 'TestTrain'
colnames(GCbyICstd)[5] = 'Set'
GCbyICstd$Cluster = factor(GCbyICstd$Cluster, levels = 0:17)
mycolors = colorRampPalette(brewer.pal(8, "Set1"))(22)

#Plot clustered transcription factor motifs
GCbyICstd_clust_plot = ggplot(data = GCbyICstd, 
                               mapping = aes(x = zstdIC, 
                                             y = zstdGC, 
                                             color = Cluster, label = TF)) + geom_point(size = 4) +
  theme_classic() +
  xlab('Information Content (z-score standardized)') +
  ylab('Motif GC% (z-score standardized)') +
  theme(axis.title.x = element_text(family = 'sans', face = 'bold', size = 22)) + 
  theme(axis.title.y = element_text(family = 'sans', face = 'bold', size = 22)) +
  theme(axis.text=element_text(size=18, color = 'black')) +
  theme(legend.text=element_text(size=18, color = 'black'), legend.title=element_text(size=18)) +
  geom_point(shape = 1,size = 4,colour = "black") + 
  scale_color_manual(values = mycolors) + 
  new_scale_color() + geom_point(aes(x = zstdIC, y = zstdGC, color = Set), alpha = 0) +
  scale_color_manual(values = c('red', 'blue')) +
  geom_text_repel(aes(label = TF),
                  size = 6, max.overlaps = 11, color = GCbyICstd$TestTrain,
                  force = 10, force_pull = 0.75, seed = 43L, min.segment.length = 0) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + coord_fixed(ratio = 1.086468)

pdf("ICbyMotifGC_TestTrain_clust.pdf", width=12, height=12)
GCbyICstd_clust_plot
dev.off()

#Plot each motif's seqlogo:
for (i in 1:length(Memes)) {
  plot.seqlogo.meme(Memes[i],
        outfile = paste(substr(Memes[i],1,nchar(Memes[i])-5),'_seqLogo.pdf',sep = ''))
}
