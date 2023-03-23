library(gridExtra)
library(grid)
#Get figure 7D (and supp) data:
system('wget https://github.com/guertinlab/Tn5bias/raw/master/Manuscript_Figures/Figure7/Figure7D_Tn5_Violin_log2Baseline.Rdata')
load('Figure7D_Tn5_Violin_log2Baseline.Rdata')

avg_baselines$baseline_log = abs(avg_baselines$baseline_log)

#All together stats mean
unscaled_mean = mean(avg_baselines[which(avg_baselines$group == 'Unscaled'),1])
seqOutBias_mean = mean(avg_baselines[which(avg_baselines$group == 'seqOutBias'),1])
RuleEnsemble_mean = mean(avg_baselines[which(avg_baselines$group == 'Rule Ensemble'),1])
paste(format(unscaled_mean, digits = 2), format(seqOutBias_mean, digits = 2), format(RuleEnsemble_mean, digits = 2))

#All together stats var
unscaled_var = var(avg_baselines[which(avg_baselines$group == 'Unscaled'),1])
seqOutBias_var = var(avg_baselines[which(avg_baselines$group == 'seqOutBias'),1])
RuleEnsemble_var = var(avg_baselines[which(avg_baselines$group == 'Rule Ensemble'),1])
paste(format(unscaled_var, digits = 2), format(seqOutBias_var, digits = 2), format(RuleEnsemble_var, digits = 2))


#All together data
unscaled_data = avg_baselines[which(avg_baselines$group == 'Unscaled'),]
seqOutBias_data = avg_baselines[which(avg_baselines$group == 'seqOutBias'),]
RuleEnsemble_data = avg_baselines[which(avg_baselines$group == 'Rule Ensemble'),]

#QQ plot unscaled
qqnorm(avg_baselines[which(avg_baselines$group == 'Unscaled'),1], pch = 1, frame = FALSE)
qqline(avg_baselines[which(avg_baselines$group == 'Unscaled'),1], col = "red", lwd = 2)

#QQ plot seqOutBias
qqnorm(avg_baselines[which(avg_baselines$group == 'seqOutBias'),1], pch = 1, frame = FALSE)
qqline(avg_baselines[which(avg_baselines$group == 'seqOutBias'),1], col = "red", lwd = 2)

#QQ plot Rule Ensemble
qqnorm(avg_baselines[which(avg_baselines$group == 'Rule Ensemble'),1], pch = 1, frame = FALSE)
qqline(avg_baselines[which(avg_baselines$group == 'Rule Ensemble'),1], col = "red", lwd = 2)

unscaled_seqOutBias_t = t.test(unscaled_data$baseline_log, seqOutBias_data$baseline_log, paired = TRUE)
unscaled_RuleEnsemble_t = t.test(unscaled_data$baseline_log, RuleEnsemble_data$baseline_log, paired = TRUE)
seqOutBias_RuleEnsemble_t = t.test(seqOutBias_data$baseline_log, RuleEnsemble_data$baseline_log, paired = TRUE)

unscaled_seqOutBias_t
unscaled_RuleEnsemble_t
seqOutBias_RuleEnsemble_t

unscaled_seqOutBias_f = var.test(unscaled_data$baseline_log, seqOutBias_data$baseline_log)
unscaled_RuleEnsemble_f = var.test(unscaled_data$baseline_log, RuleEnsemble_data$baseline_log)
seqOutBias_RuleEnsemble_f = var.test(seqOutBias_data$baseline_log, RuleEnsemble_data$baseline_log)

summary_table = data.frame(c('Unscaled', 'seqOutBias', 'Rule Ensemble'))
summary_table[,2] = c(format(unscaled_mean, digits = 2), format(seqOutBias_mean, digits = 2), format(RuleEnsemble_mean, digits = 2))
summary_table[,3] = c(format(unscaled_var, digits = 2), format(seqOutBias_var, digits = 2), format(RuleEnsemble_var, digits = 2))
summary_table[,4] = c('-', '***', '***')
summary_table[,5] = c('-', '-', '***')
summary_table[,6] = c('-', format(unscaled_seqOutBias_f[[3]], digits = 2), '**')
summary_table[,7] = c('-', '-', format(seqOutBias_RuleEnsemble_f[[3]], digits = 2))

colnames(summary_table) = c('group', 'Abs Mean', 'Abs Variance', 'Unscaled t-test p-value',
                            'seqOutBias t-test p-value', 'Unscaled F-test p-value', 'seqOutBias F-test p-value')

pdf(file = "Figure7D_summary_stats.pdf", height = 2.0, width = 19)
grid.table(summary_table, rows = rep('', nrow(summary_table)),theme=ttheme_default(base_size = 16))
grid.text("Figure 7D summary statistics", x = 0.19, y = 0.9, gp = gpar(fontsize = 20, fontface = 'bold'))
dev.off()


png(file = "Figure7D_summary_stats.png", height = 150, width = 950)
grid.table(summary_table, rows = rep('', nrow(summary_table)))
grid.text("Figure 7D summary statistics", x = 0.175, y = 0.85, gp = gpar(fontsize = 16, fontface = 'bold'))
dev.off()
