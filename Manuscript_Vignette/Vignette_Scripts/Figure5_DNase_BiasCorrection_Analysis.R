options(scipen = 100)
library(data.table)
library(lattice)
library(ggplot2)
library(grid)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
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

########################### Load in TF motif coordinates and each TF's length
Motifs = list.files('./')
Motifs = Motifs[grep('_fimo', Motifs)]

Motiflist <- vector('list', length(Motifs))
for (i in 1:length(Motifs)) {
  Motiflist[[i]] <- FIMO.to.BED(paste('./', Motifs[i], sep = ''))
}
names(Motiflist) = Motifs
#Combine minus and plus motifs
minus_Motifs = Motiflist[seq(1, 88, 2)]
plus_Motifs = Motiflist[-c(seq(1, 88, 2))]
Motiflist = vector('list', length(plus_Motifs))
for (i in 1:length(plus_Motifs)) {
  Motiflist[[i]] = rbind(plus_Motifs[[i]], minus_Motifs[[i]])
}
names(Motiflist) <- substr(names(plus_Motifs), 1, nchar(names(plus_Motifs))-22)
rm(plus_Motifs, minus_Motifs)
Motifs <- names(Motiflist)
#Get lengths for each motif
Motiflen_files = list.files('./')
Motiflen_files = Motiflen_files[grep('_fimo', Motiflen_files)]
Motiflen_files = Motiflen_files[grep('plus', Motiflen_files)]
Motiflen = NULL
Motif_len = NULL
for (i in 1:length(Motiflen_files)) {
  Motif_len = fread(Motiflen_files[i])
  Motiflen[i] = (Motif_len$stop[1] - Motif_len$start[1])+1
}
names(Motiflen) = Motifs
save(Motiflen, file = 'Motiflen.Rdata')
#Plot RE model corrected TFs
BWs <- c('DNase_Rule_Ensemble_Scaled.bigWig')

compositelist <- vector('list', length(Motifs))
for (i in 1:length(Motiflist)) {
  compositelist[[i]] <- BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                         upstream = 100, downstream = 100, 
                                         factor = Motifs[i], group = 'Rule Ensemble', ATAC = FALSE)
}

names(compositelist) <- Motifs
DNase_RE_composite = compositelist
save(DNase_RE_composite, file = 'DNase_RE_composite.Rdata')

#Plot masked TFs for comparison (this seqOutBias command was run for figure 1C)
BWs <- c('../Figure1/DNase_XXXXXXXXXXXXXXXNNNCNNXXXXXXXXXXXXXXXX.bigWig')

compositelist <- vector('list', length(Motifs))
for (i in 1:length(Motiflist)) {
  compositelist[[i]] <- BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                         upstream = 100, downstream = 100,
                                         factor = Motifs[i], group = 'seqOutBias', ATAC = FALSE)
}

names(compositelist) <- Motifs
DNase_NNNCNN_compositelist = compositelist

save(DNase_NNNCNN_compositelist, file = 'DNase_NNNCNN_compositelist.Rdata')

#Plot masked TFs for comparison (this seqOutBias command was run for figure 1C)
BWs <- c('../Figure1/DNase_Naked_unscaled.bigWig')

compositelist <- vector('list', length(Motifs))
for (i in 1:length(Motiflist)) {
  compositelist[[i]] <- BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                         upstream = 100, downstream = 100,
                                         factor = Motifs[i], group = 'Unscaled', ATAC = FALSE)
}

names(compositelist) <- Motifs
DNase_unscaled_compositelist = compositelist

save(DNase_unscaled_compositelist, file = 'DNase_unscaled_compositelist.Rdata')

##Combine rule ensemble, seqOutBias and unscaled compositelists for plotting
combined_compositelist <- vector(mode = 'list', length = 1)
for (i in 1:length(DNase_unscaled_compositelist)) {
  combined_compositelist[[i]] <- rbind(DNase_RE_composite[[i]],
                                       DNase_unscaled_compositelist[[i]],
                                       DNase_NNNCNN_compositelist[[i]])
}
names(combined_compositelist) <- names(DNase_unscaled_compositelist)

#Separate out the TFs to be plotted (for 5C) and their motif lengths
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



###Subset test motifs out
test_DNase_RE_compositelist = DNase_RE_composite[which(names(DNase_RE_composite) %in% test)]
test_unscaled_compositelist = 
  DNase_unscaled_compositelist[which(names(DNase_unscaled_compositelist) %in% test)]
test_DNase_NNNCNN_compositelist = 
  DNase_NNNCNN_compositelist[which(names(DNase_NNNCNN_compositelist) %in% test)]
Motiflen = Motiflen[which(names(Motiflen) %in% test)]
#Make data frame for single nucleotide measurements of each treatment
test_DNase_RE = do.call(rbind, test_DNase_RE_compositelist)
singlenuc_frame = NULL
singlenuc_store = NULL
for (i in 1:length(Motiflen)) {
  singlenuc_store = test_DNase_RE[which(test_DNase_RE$factor == names(Motiflen)[i] &
                                          test_DNase_RE$x < ((Motiflen[i]/2)+10) &
                                          test_DNase_RE$x > - ((Motiflen[i]/2)+10)),]
  singlenuc_frame = rbind(singlenuc_frame, singlenuc_store)
}
test_DNase_RE = singlenuc_frame

test_unscaled = do.call(rbind, test_unscaled_compositelist)
singlenuc_frame = NULL
singlenuc_store = NULL
for (i in 1:length(Motiflen)) {
  singlenuc_store = test_unscaled[which(test_unscaled$factor == names(Motiflen)[i] &
                                          test_unscaled$x < ((Motiflen[i]/2)+10) &
                                          test_unscaled$x > - ((Motiflen[i]/2)+10)),]
  singlenuc_frame = rbind(singlenuc_frame, singlenuc_store)
}
test_unscaled = singlenuc_frame

test_DNase_NNNCNN = do.call(rbind, test_DNase_NNNCNN_compositelist)
singlenuc_frame = NULL
singlenuc_store = NULL
for (i in 1:length(Motiflen)) {
  singlenuc_store = test_DNase_NNNCNN[which(test_DNase_NNNCNN$factor == names(Motiflen)[i] &
                                              test_DNase_NNNCNN$x < ((Motiflen[i]/2)+10) &
                                              test_DNase_NNNCNN$x > - ((Motiflen[i]/2)+10)),]
  singlenuc_frame = rbind(singlenuc_frame, singlenuc_store)
}
test_DNase_NNNCNN = singlenuc_frame

#to get log2 diff- first change to percent diff from baseline
#Calculate distance between target and output
DNase_RE_log = NULL
DNase_NNNCNN_log = NULL
unscaled_log = NULL
for (i in 1:length(test_unscaled$factor)) {
  DNase_RE_log[i] = log2(test_DNase_RE$est[i] / 0.006831447)
  unscaled_log[i] = log2(test_unscaled$est[i] / 0.006831447)
  DNase_NNNCNN_log[i] = log2(test_DNase_NNNCNN$est[i] / 0.006831447)
}

unscaled_log = data.frame(unscaled_log)
unscaled_log[,2] = test_unscaled$factor
unscaled_log[,3] = 'Unscaled'
colnames(unscaled_log) = c('Difference', 'Factor', 'Treatment')

DNase_RE_log = data.frame(DNase_RE_log)
DNase_RE_log[,2] = test_unscaled$factor
DNase_RE_log[,3] = 'Rule Ensemble'
colnames(DNase_RE_log) = c('Difference', 'Factor', 'Treatment')

DNase_NNNCNN_log = data.frame(DNase_NNNCNN_log)
DNase_NNNCNN_log[,2] = test_unscaled$factor
DNase_NNNCNN_log[,3] = 'seqOutBias'
colnames(DNase_NNNCNN_log) = c('Difference', 'Factor', 'Treatment')

log_fold = rbind(DNase_NNNCNN_log, DNase_RE_log, unscaled_log)
log_fold$Treatment = factor(log_fold$Treatment, levels = c('Unscaled',
                                                           'seqOutBias',
                                                           'Rule Ensemble'))

x = log_fold

x[,3] <- as.factor(x[,3])
x[,3] = factor(x[,3], levels = c('Unscaled', 'seqOutBias', 'Rule Ensemble'))

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
dot_frame = test_unscaled
dot_frame = dot_frame[,-c(5)]
colnames(dot_frame) = c('Unscaled_x', 'Unscaled_est', 'Unscaled_factor', 'Unscaled_group')

dot_frame = cbind(dot_frame, test_DNase_RE[,1:4])
colnames(dot_frame)[5:8] = c('RE_x', 'RE_est', 'RE_factor', 'RE_group')
dot_frame = cbind(dot_frame, test_DNase_NNNCNN[,1:4])
colnames(dot_frame)[9:12] = c('sOB_x', 'sOB_est', 'sOB_factor', 'sOB_group')
identical(dot_frame$Unscaled_factor, dot_frame$RE_factor)

dot_frame = dot_frame[,-c(4,5,7,8,9,11,12)]
dot_frame$calc_avg = 0.006831447
dot_frame$RE_absdiff = abs(dot_frame$RE_est - dot_frame$calc_avg)
dot_frame$sOB_absdiff = abs(dot_frame$sOB_est - dot_frame$calc_avg)
rownames(dot_frame) = 1:nrow(dot_frame)

dot_frame$sOB_sub_RE = dot_frame$sOB_absdiff - dot_frame$RE_absdiff
dot_frame$plot_col  = 'RE Improvement over sOB'

length(which(dot_frame$sOB_absdiff > dot_frame$RE_absdiff))/ nrow(dot_frame)

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
