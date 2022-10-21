library(lattice)
library(data.table)
library(grid)
library(ggplot2)
setwd('../output')
################################################################################
################################################################################
################################################################################
load('../data/220904_Tn5_singlenuc_composites.Rdata')
################################################################################

dot_frame = singlenuc_frame[which(singlenuc_frame$group == 'Unscaled'),]
dot_frame = dot_frame[,-c(5)]
colnames(dot_frame) = c('Unscaled_x', 'Unscaled_est', 'Unscaled_factor', 'Unscaled_group')

dot_frame = cbind(dot_frame, singlenuc_frame[which(singlenuc_frame$group == 'Rule Ensemble'),1:4])
colnames(dot_frame)[5:8] = c('RE_x', 'RE_est', 'RE_factor', 'RE_group')
dot_frame = cbind(dot_frame, singlenuc_frame[which(singlenuc_frame$group == 'seqOutBias'),1:4])
colnames(dot_frame)[9:12] = c('sOB_x', 'sOB_est', 'sOB_factor', 'sOB_group')
identical(dot_frame$Unscaled_x, dot_frame$RE_x)


dot_frame = dot_frame[,-c(4,5,7,8,9,11,12)]
dot_frame$calc_avg = 0.006140407
dot_frame$RE_absdiff = abs(dot_frame$RE_est - dot_frame$calc_avg)
dot_frame$sOB_absdiff = abs(dot_frame$sOB_est - dot_frame$calc_avg)
rownames(dot_frame) = 1:nrow(dot_frame)





plot(dot_frame$Unscaled_x[1:35], dot_frame$Unscaled_est[1:35], type = 'l')
lines(dot_frame$Unscaled_x[1:35], dot_frame$sOB_est[1:35], col = 'red')
lines(dot_frame$Unscaled_x[1:35], dot_frame$RE_est[1:35], col = 'blue')


dot_frame$sOB_sub_RE = dot_frame$sOB_absdiff - dot_frame$RE_absdiff
dot_frame$plot_col  = 'RE Improvement over sOB'


ggplot(dot_frame, aes(x=sOB_sub_RE, y = plot_col)) + 
  geom_density(position = 'identity')


d <- density(dot_frame$sOB_sub_RE)
plot(d)

pdf('Figure6E_DNase_sOB_RE_improvement_BW.pdf', width=10, height=10)
ggplot(dot_frame, aes(x=plot_col, y = sOB_sub_RE)) + 
  geom_violin(trim = FALSE, color = 'black', fill = 'light blue') +  xlab("") + ylab("Rule Ensemble Improvement") + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.4) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text = element_text(size = 36), axis.title = element_text(size = 42, family = 'sans', face = 'bold'),
        axis.text.y = element_text(colour = "black", family = 'sans')) + geom_abline(slope = 0, lwd = 2, color = 'red', linetype = "dashed")
dev.off()

