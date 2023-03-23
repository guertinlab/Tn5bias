library(gridExtra)
library(grid)
#Get figure 6C (and supp) data:
system('wget https://github.com/guertinlab/Tn5Bias/raw/master/Manuscript_Figures/Figure6/Figure6C_Tn5_singlenuc_log2_comparison.Rdata')
system('wget https://github.com/guertinlab/Tn5bias/raw/master/Manuscript_Figures/Figure6/Supplemental_Tn5_singlenuc_log2_comparison.Rdata')
load('Figure6C_Tn5_singlenuc_log2_comparison.Rdata')
load('Supplemental_Tn5_singlenuc_log2_comparison.Rdata')

singlenuc_log = rbind(singlenuc_log, supp_singlenuc_log)
singlenuc_log$Difference = abs(singlenuc_log$Difference)

RE_frame = data.frame()
seqOutBias_frame = data.frame()
unscaled_frame = data.frame()

for (i in 1:length(unique(singlenuc_log$Factor))) {
  RE_frame[i,1] = mean(singlenuc_log[which(singlenuc_log$Factor == unique(singlenuc_log$Factor)[i] & singlenuc_log$Treatment == 'Rule Ensemble'),1])
  RE_frame[i,2] = var(singlenuc_log[which(singlenuc_log$Factor == unique(singlenuc_log$Factor)[i] & singlenuc_log$Treatment == 'Rule Ensemble'),1])
  RE_frame[i,3] = unique(singlenuc_log$Factor)[i]
  seqOutBias_frame[i,1] = mean(singlenuc_log[which(singlenuc_log$Factor == unique(singlenuc_log$Factor)[i] & singlenuc_log$Treatment == 'seqOutBias'),1])
  seqOutBias_frame[i,2] = var(singlenuc_log[which(singlenuc_log$Factor == unique(singlenuc_log$Factor)[i] & singlenuc_log$Treatment == 'seqOutBias'),1])
  seqOutBias_frame[i,3] = unique(singlenuc_log$Factor)[i]
  unscaled_frame[i,1] = mean(singlenuc_log[which(singlenuc_log$Factor == unique(singlenuc_log$Factor)[i] & singlenuc_log$Treatment == 'Unscaled'),1])
  unscaled_frame[i,2] = var(singlenuc_log[which(singlenuc_log$Factor == unique(singlenuc_log$Factor)[i] & singlenuc_log$Treatment == 'Unscaled'),1])
  unscaled_frame[i,3] = unique(singlenuc_log$Factor)[i]
}

colnames(unscaled_frame) = c('mean', 'var', 'TF')
colnames(seqOutBias_frame) = c('mean', 'var', 'TF')
colnames(RE_frame) = c('mean', 'var', 'TF')

combined_frame = data.frame(unscaled_frame$TF)
combined_frame[,2:7] =  c(unscaled_frame$mean, seqOutBias_frame$mean, RE_frame$mean,
                       unscaled_frame$var, seqOutBias_frame$var, RE_frame$var)
combined_frame = format(combined_frame, digits = 2)

colnames(combined_frame) = c('Transcription Factor', 'Unscaled Abs Mean', 'seqOutBias Abs Mean', 'Rule Ensemble Abs Mean',
                             'Unscaled Abs Variance', 'seqOutBias Abs Variance', 'Rule Ensemble Abs Variance')


pdf(file = "Figure6C_TF_stats.pdf", height = 6, width = 15)
grid.table(combined_frame, rows = rep('', nrow(combined_frame)))
dev.off()

#All together stats mean
unscaled_mean = mean(singlenuc_log[which(singlenuc_log$Treatment == 'Unscaled'),1])
seqOutBias_mean = mean(singlenuc_log[which(singlenuc_log$Treatment == 'seqOutBias'),1])
RuleEnsemble_mean = mean(singlenuc_log[which(singlenuc_log$Treatment == 'Rule Ensemble'),1])
paste(format(unscaled_mean, digits = 2), format(seqOutBias_mean, digits = 2), format(RuleEnsemble_mean, digits = 2))

#All together stats var
unscaled_var = var(singlenuc_log[which(singlenuc_log$Treatment == 'Unscaled'),1])
seqOutBias_var = var(singlenuc_log[which(singlenuc_log$Treatment == 'seqOutBias'),1])
RuleEnsemble_var = var(singlenuc_log[which(singlenuc_log$Treatment == 'Rule Ensemble'),1])
paste(format(unscaled_var, digits = 2), format(seqOutBias_var, digits = 2), format(RuleEnsemble_var, digits = 2))

#QQ plot unscaled
qqnorm(singlenuc_log[which(singlenuc_log$Treatment == 'Unscaled'),1], pch = 1, frame = FALSE)
qqline(singlenuc_log[which(singlenuc_log$Treatment == 'Unscaled'),1], col = "red", lwd = 2)

#QQ plot seqOutBias
qqnorm(singlenuc_log[which(singlenuc_log$Treatment == 'seqOutBias'),1], pch = 1, frame = FALSE)
qqline(singlenuc_log[which(singlenuc_log$Treatment == 'seqOutBias'),1], col = "red", lwd = 2)

#QQ plot Rule Ensemble
qqnorm(singlenuc_log[which(singlenuc_log$Treatment == 'Rule Ensemble'),1], pch = 1, frame = FALSE)
qqline(singlenuc_log[which(singlenuc_log$Treatment == 'Rule Ensemble'),1], col = "red", lwd = 2)

#All together data
unscaled_data = singlenuc_log[which(singlenuc_log$Treatment == 'Unscaled'),]
seqOutBias_data = singlenuc_log[which(singlenuc_log$Treatment == 'seqOutBias'),]
RuleEnsemble_data = singlenuc_log[which(singlenuc_log$Treatment == 'Rule Ensemble'),]

unscaled_seqOutBias_w = wilcox.test(unscaled_data$Difference, seqOutBias_data$Difference, paired = TRUE)
unscaled_RuleEnsemble_w = wilcox.test(unscaled_data$Difference, RuleEnsemble_data$Difference, paired = TRUE)
seqOutBias_RuleEnsemble_w = wilcox.test(seqOutBias_data$Difference, RuleEnsemble_data$Difference, paired = TRUE)

unscaled_seqOutBias_w
unscaled_RuleEnsemble_w
seqOutBias_RuleEnsemble_w

unscaled_seqOutBias_f = var.test(unscaled_data$Difference, seqOutBias_data$Difference)
unscaled_RuleEnsemble_f = var.test(unscaled_data$Difference, RuleEnsemble_data$Difference)
seqOutBias_RuleEnsemble_f = var.test(seqOutBias_data$Difference, RuleEnsemble_data$Difference)

summary_table = data.frame(c('Unscaled', 'seqOutBias', 'Rule Ensemble'))
summary_table[,2] = c(format(unscaled_mean, digits = 2), format(seqOutBias_mean, digits = 2), format(RuleEnsemble_mean, digits = 2))
summary_table[,3] = c(format(unscaled_var, digits = 2), format(seqOutBias_var, digits = 2), format(RuleEnsemble_var, digits = 2))
summary_table[,4] = c('-', '***', '***')
summary_table[,5] = c('-', '-', '***')
summary_table[,6] = c('-', '***', '***')
summary_table[,7] = c('-', '-', '***')

colnames(summary_table) = c('Treatment', 'Abs Mean', 'Abs Variance', 'Unscaled Mann-Whitney U test p-value',
                            'seqOutBias Mann-Whitney U test p-value', 'Unscaled F-test p-value', 'seqOutBias F-test p-value')

pdf(file = "Figure6C_summary_stats.pdf", height = 2.0, width = 20)
grid.table(summary_table, rows = rep('', nrow(summary_table)),theme=ttheme_default(base_size = 16))
grid.text("Figure 6C summary statistics", x = 0.122, y = 0.9, gp = gpar(fontsize = 20, fontface = 'bold'))
dev.off()


png(file = "Figure6C_summary_stats.png", height = 150, width = 1150)
grid.table(summary_table, rows = rep('', nrow(summary_table)))
grid.text("Figure 6C summary statistics", x = 0.155, y = 0.85, gp = gpar(fontsize = 16, fontface = 'bold'))
dev.off()
