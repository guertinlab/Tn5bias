library(ggplot2)
load('Figure6D_RuleEnsemble_Improvement.Rdata')
################################################################################

pdf('Figure6D_Tn5_sOB_RE_improvement_BW.pdf', width=10, height=12)
ggplot(dot_frame, aes(x=plot_col, y = sOB_sub_RE)) + 
  geom_violin(trim = FALSE, color = 'black', fill = 'light blue') +
  xlab("") + ylab("Rule Ensemble Improvement") + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.4) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text = element_text(size = 36), axis.title = element_text(size = 42,
        family = 'sans', face = 'bold'),
        axis.text.y = element_text(colour = "black", family = 'sans')) +
        geom_abline(slope = 0, lwd = 2, color = 'red', linetype = "dashed")
dev.off()
