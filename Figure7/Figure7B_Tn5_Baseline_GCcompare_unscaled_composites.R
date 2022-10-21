library(ggplot2)
library(ggrepel)
setwd('../output')
load('../data/Figure7B_Tn5_GCcon_baseline.Rdata')

pdf("Figure7B_Tn5_unscaled_BaselineVsGCcon.pdf", width=6, height=6)
ggplot(Tn5_GCcon_baseline, aes(x = V1, y = V3, label = V4)) +
  geom_point(stat = 'identity') +
  xlim(0, 0.008) +
  ylim(0, 0.8) +
  theme_classic() + xlab('Baseline Signal') + ylab('Motif GC%') + 
  #geom_text_repel(size = 2, stat = 'identity', force_pull = 0.75, force = 5, max.overlaps = 20, seed = 42) + 
  theme(axis.text=element_text(size=22, color = 'black'),
        axis.title=element_text(size=30, color = 'black'),
        axis.text.y=element_text(angle = 90, hjust = 0.4),
        axis.text.x=element_text(hjust = 0.60)) 
dev.off()
