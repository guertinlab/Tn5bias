library(ggplot2)
library(lattice)
setwd('../output')
library(data.table)
library(ggplot2)
library(scclust)
library(RColorBrewer)
library(ggrepel)
library(ggnewscale)
################################################################################
################################################################################
#GCbyICstd$Cluster = factor(GCbyICstd$Cluster, levels = seq(0, max(GCbyICstd$Cluster),1))
################################################################################
mycolors = colorRampPalette(brewer.pal(8, "Set1"))(22)
load('../data/Figure3A_ICbyMotifGC_TestTrain_clust.Rdata')
GCbyICstd_kmeans_plot = ggplot(data = GCbyICstd, 
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
                  size = 6, max.overlaps = 11, color = GCbyICstd$TestTrain,  force = 10, force_pull = 0.75, seed = 43L, min.segment.length = 0) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + coord_fixed(ratio = 1.086468)

pdf("ICbyMotifGC_TestTrain_clust.pdf", width=12, height=12)
GCbyICstd_kmeans_plot
dev.off()
