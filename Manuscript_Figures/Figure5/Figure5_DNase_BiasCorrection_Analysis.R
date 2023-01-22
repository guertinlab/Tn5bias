library(grid)
library(lattice)
library(ggplot2)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
#Set test and training sets
train = c( 'ASCL1', 'ATF3', 'CLOCK', 'CTCF',
           'DLX2', 'EBF1', 'EGR1', 'ELF1',
           'FOS', 'FOXA1', 'FOXK2', 'GATA2',
           'GLIS1', 'IRF1', 'JUN', 'LEF1',
           'MEIS2', 'MLX', 'MYC', 'PPARG', 
           'RUNX1', 'SPIB', 'SRY', 'STAT1', 'TGIF1')
test = c('AR', 'CEBPB', 'E2F1', 'ESR1',
         'FERD3L', 'HOXC12', 'HSF1', 'MAX',
         'MEF2A', 'MGA', 'NFATC3', 'NR2F2',
         'POU3F1', 'REST', 'Six3', 'SP1',
         'USF1', 'TEAD1')

load('Figure5C_composites.Rdata')
load('Figure5C_Motiflen.Rdata')
#Separate out the TFs to be plotted (for 3C) and their motif lengths
figure5C_plot = do.call(rbind, combined_compositelist[c(4, 20, 41)])

mlen <- Motiflen[c(4, 20, 41)]
mlen <- as.data.frame(mlen)
rownames(mlen) <- levels(as.factor(figure5C_plot$factor))


plot.composites(figure5C_plot, legend = TRUE, 
                pdf_name = 'Figure5C_DNasePN_PRE_NNNCNN_comparison',
                figwidth = 8, figheight = 4,
                ylabel = '',
                xlabel = '',
                motifline = TRUE, Motiflen = mlen, layoutgrid = c(3,1),
                x_axis_range = -20:20, X_axis_ticks = seq(-20,20,10),
                Y_axis_ticks = seq(0,0.05,0.01), nYaxisdigits = 2,
                hline_val = 0.006831447, y_axis_range = seq(0,0.02,0.01),
                y_axis = FALSE, hline = TRUE, Y_ticks = FALSE, labsize = 0.85)

#Plot supplemental composites
supp_figure5C_plot = do.call(rbind, combined_compositelist[-c(4, 20, 41)])
supp_figure5C_plot = supp_figure5C_plot[which(supp_figure5C_plot$factor %in% test),]

mlen = Motiflen[-c(4, 20, 41)]
mlen = mlen[which(names(mlen) %in% test)]
mlen = as.data.frame(mlen)
rownames(mlen) = levels(as.factor(supp_figure5C_plot$factor))


plot.composites(supp_figure5C_plot, legend = TRUE, 
                pdf_name = 'Supplemental5B_DNasePN_PRE_NNNCNN_comparison',
                figwidth = 12, figheight = 10,
                ylabel = '',
                xlabel = '',
                motifline = TRUE, Motiflen = mlen, layoutgrid = c(5,3),
                x_axis_range = -20:20, X_axis_ticks = seq(-20,20,10),
                Y_ticks = FALSE,
                hline_val = 0.006831447, y_axis_range = seq(0,0.02,0.01),
                y_axis = FALSE, hline = TRUE, labsize = 0.85)



###Log2 comparison
load('Figure5B_log2_comparison.Rdata')

#Welch t-test to see if treatments' outputs are different
t.test(x[which(x$Treatment=='seqOutBias'),1], x[which(x$Treatment=='Rule Ensemble'),1])
t.test(x[which(x$Treatment=='Unscaled'),1], x[which(x$Treatment=='Rule Ensemble'),1])
t.test(x[which(x$Treatment=='Unscaled'),1], x[which(x$Treatment=='seqOutBias'),1])

#F-test to see if treatments' variances are different for each TF
vartests = data.frame()
for (i in 1:length(unique(x$Factor))) {
  vartests[i,1] = as.numeric(var.test(x[which(x$Treatment=='seqOutBias' & x$Factor==unique(x$Factor)[i]),1],
                         x[which(x$Treatment=='Unscaled' & x$Factor==unique(x$Factor)[i]),1])[3])
  vartests[i,2] = as.numeric(var.test(x[which(x$Treatment=='Rule Ensemble' & x$Factor==unique(x$Factor)[i]),1],
                         x[which(x$Treatment=='Unscaled' & x$Factor==unique(x$Factor)[i]),1])[3])
  vartests[i,3] = as.numeric(var.test(x[which(x$Treatment=='Rule Ensemble' & x$Factor==unique(x$Factor)[i]),1],
                         x[which(x$Treatment=='seqOutBias' & x$Factor==unique(x$Factor)[i]),1])[3])
}
colnames(vartests) = c('seqOutBias Unscaled', 'Rule Ensemble Unscaled', 'Rule Ensemble seqOutBias')
rownames(vartests)  = unique(x$Factor)


#Subset x to only include top 9
supp_x = x[-c(which(x$Factor == 'CEBPB' | x$Factor == 'E2F1' | x$Factor == 'HOXC12' |
                      x$Factor == 'HSF1' | x$Factor == 'MEF2A' | x$Factor == 'NFATC3' |
                      x$Factor == 'POU3F1' | x$Factor == 'Six3' | x$Factor == 'TEAD1')),]
x = x[c(which(x$Factor == 'CEBPB' | x$Factor == 'E2F1' | x$Factor == 'HOXC12' |
                x$Factor == 'HSF1' | x$Factor == 'MEF2A' | x$Factor == 'NFATC3' |
                x$Factor == 'POU3F1' | x$Factor == 'Six3' | x$Factor == 'TEAD1')),]

log2plot = bwplot(Difference ~ Treatment | Factor , data = x,
                  between=list(y=0.5, x = 0.5),
                  scales=list(x=list(draw=TRUE),rot = 45, alternating=c(1,1,1,1),cex=1.2,font=1),
                  #xlab = '',
                  ylim = c(-2.10, 2.10),
                  ylab =expression("log"[2]~frac("Unbiased Signal","Model Signal Output")),
                  horizontal =FALSE,  col= 'black',
                  aspect = 2,
                  par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                                    par.ylab.text=list(cex=1.2,font=1),
                                    par.main.text=list(cex=1.2, font=1),
                                    plot.symbol = list(col='black', lwd=0, pch =19, cex = 0.4)),
                  par.strip.text = list(cex = 1.2),
                  strip = function(..., which.panel, bg) {
                    bg.col = c("white")
                    strip.default(..., which.panel = which.panel,
                                  bg = rep(bg.col, length = which.panel)[which.panel])
                  },
                  panel = function(..., box.ratio, col) {
                    panel.abline(h = 0, col = 'grey45', lty = 2)
                    panel.violin.hack(..., col = c("#00000090", "#FF0000", "#0000FF"),
                                      varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                    panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                    panel.bwplot(..., pch = '|', do.out = FALSE)
                    
                  })

pdf('Figure5B_DNase_log2_comparison.pdf', useDingbats = FALSE, width=10.83, height=6)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939300", lwd=2),
                box.rectangle = list(col = '#93939300', lwd=1.6),
                plot.symbol = list(col='#93939300', lwd=1.6, pch ='.'))

################################################################################
print(update(log2plot, layout=c(9,1)))
dev.off()

log2plot = bwplot(Difference ~ Treatment | Factor , data = supp_x,
                  between=list(y=0.5, x = 0.5),
                  scales=list(x=list(draw=TRUE),rot = 45, alternating=c(1,1,1,1),cex=1.2,font=1),
                  #xlab = '',
                  ylim = c(-2.10, 2.10),
                  ylab =expression("log"[2]~frac("Unbiased Signal","Model Signal Output")),
                  horizontal =FALSE,  col= 'black',
                  aspect = 2,
                  par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                                    par.ylab.text=list(cex=1.2,font=1),
                                    par.main.text=list(cex=1.2, font=1),
                                    plot.symbol = list(col='black', lwd=0, pch =19, cex = 0.4)),
                  par.strip.text = list(cex = 1.2),
                  strip = function(..., which.panel, bg) {
                    bg.col = c("white")
                    strip.default(..., which.panel = which.panel,
                                  bg = rep(bg.col, length = which.panel)[which.panel])
                  },
                  panel = function(..., box.ratio, col) {
                    panel.abline(h = 0, col = 'grey45', lty = 2)
                    panel.violin.hack(..., col = c("#00000090", "#FF0000", "#0000FF"),
                                      varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                    panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                    panel.bwplot(..., pch = '|', do.out = FALSE)
                    
                  })

pdf('Supplemental_Figure5B_DNase_log2_comparison.pdf', useDingbats = FALSE, width=10.83, height=6)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939300", lwd=2),
                box.rectangle = list(col = '#93939300', lwd=1.6),
                plot.symbol = list(col='#93939300', lwd=1.6, pch ='.'))

################################################################################
print(update(log2plot, layout=c(9,1)))
dev.off()

#Plot Improvement of rule ensemble over seqOutBias
load('Figure5D_DNase_RuleEnsemble_improvement.Rdata')

pdf('Figure5D_DNase_sOB_RE_improvement_BW.pdf', width=10, height=12)
ggplot(dot_frame, aes(x=plot_col, y = sOB_sub_RE)) + 
  geom_violin(trim = FALSE, color = 'black', fill = 'light blue') +
  xlab("") + ylab("Rule Ensemble Improvement") + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.4) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text = element_text(size = 36),
        axis.title = element_text(size = 42, family = 'sans', face = 'bold'),
        axis.text.y = element_text(colour = "black", family = 'sans')) +
        geom_abline(slope = 0, lwd = 2, color = 'red', linetype = "dashed")
dev.off()
