source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
library(lattice)
library(grid)
################################################################################
load('../Figure6/Tn5_unscaled_compositelist.Rdata')
load('../Figure6/Tn5_RE_compositelist.Rdata')
load('../Figure6/Tn5_NNNCNNNN_compositelist.Rdata')
load('../Figure6/test_Motiflen.Rdata')
names(Motiflen) = names(Tn5_unscaled_compositelist)
#Remove single nucleotide values
unscaled_composites = rbind(do.call(rbind, Tn5_RE_compositelist), 
                            do.call(rbind, Tn5_unscaled_compositelist), 
                            do.call(rbind, Tn5_NNNCNNNN_compositelist))
unscaled_composites$group[which(unscaled_composites$group == 'NNNCNNNN')] = 'seqOutBias'

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

#Determine log2 of baseline values
baseline_log = NULL
for (i in 1:length(baseline_frame$x)) {
  baseline_log[i] = log2(baseline_frame$est[i] / 0.006140407)
}
baseline_log = data.frame(baseline_log)
baseline_log[,2:4] = baseline_frame[,c(1,3:4)]

#Determine upstream and downstream averages for each group
upstream_baseline_RE_avg = NULL
downstream_baseline_RE_avg = NULL

upstream_baseline_sOB_avg = NULL
downstream_baseline_sOB_avg = NULL

upstream_baseline_unscaled_avg = NULL
downstream_baseline_unscaled_avg = NULL

for (i in 1:length(unique(baseline_log$factor))) {
  upstream_baseline_RE_avg[i] = mean(baseline_log$baseline_log[which(baseline_log$factor ==
                                unique(baseline_log$factor)[i] & 
                                baseline_log$group == 'Rule Ensemble' & baseline_log$x < 0)])
  downstream_baseline_RE_avg[i] = mean(baseline_log$baseline_log[which(baseline_log$factor ==
                                unique(baseline_log$factor)[i] & 
                                baseline_log$group == 'Rule Ensemble' & baseline_log$x > 0)])
  
  
  upstream_baseline_sOB_avg[i] = mean(baseline_log$baseline_log[which(baseline_log$factor ==
                                unique(baseline_log$factor)[i] & 
                                baseline_log$group == 'seqOutBias' & baseline_log$x < 0)])
  downstream_baseline_sOB_avg[i] = mean(baseline_log$baseline_log[which(baseline_log$factor ==
                                unique(baseline_log$factor)[i] & 
                                baseline_log$group == 'seqOutBias' & baseline_log$x > 0)])
  
  upstream_baseline_unscaled_avg[i] = mean(baseline_log$baseline_log[which(baseline_log$factor ==
                                unique(baseline_log$factor)[i] & 
                                baseline_log$group == 'Unscaled' & baseline_log$x < 0)])
  downstream_baseline_unscaled_avg[i] = mean(baseline_log$baseline_log[which(baseline_log$factor ==
                                unique(baseline_log$factor)[i] & 
                                baseline_log$group == 'Unscaled' & baseline_log$x > 0)])
  
}
baseline_RE_avg = as.data.frame(c(upstream_baseline_RE_avg,
                                  downstream_baseline_RE_avg))
baseline_sOB_avg = as.data.frame(c(upstream_baseline_sOB_avg,
                                   downstream_baseline_sOB_avg))
baseline_unscaled_avg = as.data.frame(c(upstream_baseline_unscaled_avg,
                                        downstream_baseline_unscaled_avg))

baseline_unscaled_avg[,2] = 'Unscaled'
baseline_sOB_avg[,2] = 'seqOutBias'
baseline_RE_avg[,2] = 'Rule Ensemble'

colnames(baseline_unscaled_avg) = c('baseline_log', 'group')
colnames(baseline_sOB_avg) = c('baseline_log', 'group')
colnames(baseline_RE_avg) = c('baseline_log', 'group')

avg_baselines = rbind(baseline_unscaled_avg,baseline_sOB_avg,baseline_RE_avg)
avg_baselines$group = factor(avg_baselines$group,
                             levels = c('Unscaled', 'seqOutBias', 'Rule Ensemble'))

#plot baselines:
pdf('Figure7D_Tn5_baseline_log2_comparison.pdf', useDingbats = FALSE, width=10.83, height=6)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939300", lwd=2),
                box.rectangle = list(col = '#93939300', lwd=1.6),
                plot.symbol = list(col='#93939300', lwd=1.6, pch ='.'))

################################################################################
print(bwplot(baseline_log ~ group , data = avg_baselines,
             between=list(y=0.5, x = 0.5),
             scales=list(x=list(draw=TRUE),rot = 45, alternating=c(1,1,1,1),cex=1.2,font=1),
             #xlab = '',
             ylim = c(-2.5, 2.5),
             ylab =expression("log"[2]~frac("Unbiased Signal","Model Signal Output")),
             horizontal =FALSE,  col= 'black',
             aspect = 2,
             par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                               par.ylab.text=list(cex=1.5,font=1),
                               par.main.text=list(cex=1.2, font=1),
                               plot.symbol = list(col='black', lwd=0, pch =19, cex = 0.5)),
             par.strip.text = list(cex = 1.2),
             strip = function(..., which.panel, bg) {
               bg.col = c("grey90")
               strip.default(..., which.panel = which.panel,
                             bg = rep(bg.col, length = which.panel)[which.panel])
             },
             panel = function(..., box.ratio, col) {
               panel.abline(h = 0, col = 'grey45', lty = 2)
               panel.violin.hack(..., col = c("#A1A3AB", "#FF0000", "#0000FF"),
               varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
               panel.stripplot(..., col='#54545380', do.out=FALSE,
                               jitter.data=TRUE, amount = 0.2, pch = 16)
               panel.bwplot(..., pch = '|', do.out = FALSE)
               
             }))
dev.off()
