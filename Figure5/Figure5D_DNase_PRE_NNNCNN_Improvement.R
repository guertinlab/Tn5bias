library(ggplot2)
library(lattice)
library(grid)
setwd('../output')
################################################################################
################################################################################
################################################################################
load('../data/Figure5D_DNase_improvement.Rdata')

pdf('Figure5D_DNase_sOB_RE_improvement_BW.pdf', width=5, height=12)
ggplot(dot_frame, aes(x=plot_col, y = sOB_sub_RE)) + 
  geom_violin(trim = FALSE, color = 'black', fill = 'light blue') +  xlab("") + ylab("Rule Ensemble Improvement") + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.4) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text = element_text(size = 32), axis.title = element_text(size = 40, family = 'sans', face = 'bold'),
        axis.text.y = element_text(colour = "black", family = 'sans')) + geom_abline(slope = 0, lwd = 2, color = 'red', linetype = "dashed")
dev.off()
