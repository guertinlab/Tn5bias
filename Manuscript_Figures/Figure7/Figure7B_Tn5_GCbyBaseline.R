library(ggplot2)
load('Figure7B_BaselinevsGCcon.Rdata')
################################################################################

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

