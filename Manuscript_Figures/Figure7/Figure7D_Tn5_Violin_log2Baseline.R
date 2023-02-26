source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
library(lattice)
library(grid)
################################################################################
load('Figure7D_Tn5_Violin_log2Baseline.Rdata')
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
