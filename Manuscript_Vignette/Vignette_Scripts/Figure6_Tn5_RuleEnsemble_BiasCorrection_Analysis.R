library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
library(bigWig)
library(lattice)
library(grid)
################################################################################
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
################################################################################
#Load in genomic locations determined using FIMO
Motifs <- list.files('../Figure5')
Motifs = Motifs[grep('_fimo', Motifs)]
Motiflist <- vector('list', length(Motifs))
for (i in 1:length(Motifs)) {
  Motiflist[[i]] <- FIMO.to.BED(paste('../Figure5/',Motifs[i], sep = ''))
}
names(Motiflist) = Motifs
#Combine the minus and plus FIMO motifs
minus_Motifs = Motiflist[seq(1, 88, 2)]
plus_Motifs = Motiflist[-c(seq(1, 88, 2))]
Motiflist = vector('list', length(plus_Motifs))
for (i in 1:length(plus_Motifs)) {
  Motiflist[[i]] = rbind(plus_Motifs[[i]], minus_Motifs[[i]])
}
names(Motiflist) = substr(names(plus_Motifs), 1 , nchar(names(plus_Motifs))-22)
rm(plus_Motifs, minus_Motifs)
Motiflist = Motiflist[c(which(names(Motiflist) %in% test))]
#Make a vector of each motif's length (for graphing)
Motiflen_files = list.files('../Figure5')
Motiflen_files = Motiflen_files[grep('_fimo', Motiflen_files)]
Motiflen_files = Motiflen_files[grep('plus', Motiflen_files)]
Motiflen_files = 
  Motiflen_files[which(substr(Motiflen_files,1, nchar(Motiflen_files)-22) %in% test)]
Motiflen = NULL
Motif_len = NULL
for (i in 1:length(Motiflen_files)) {
  Motif_len = fread(paste( '../Figure5/',Motiflen_files[i], sep = ''))
  Motiflen[i] = (Motif_len$stop[1] - Motif_len$start[1])+1
}
names(Motiflen) = names(Motiflist)
save(Motiflen, file = 'test_Motiflen.Rdata')
rm(Motif_len)
#################################################################################
#Make composite plots of unscaled values 
BWs <- c('../Figure1/C1_gDNA_rep1.bigWig')
compositelist <- vector('list', length(Motiflist))
for (i in 1:length(Motiflist)) {
  compositelist[[i]] = BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                        upstream = 100, downstream = 100,
                                        factor = names(Motiflist[i]),
                                        group = 'Unscaled', ATAC = TRUE)
}
names(compositelist) = names(Motiflist)
Tn5_unscaled_compositelist = compositelist
save(Tn5_unscaled_compositelist, file = 'Tn5_unscaled_compositelist.Rdata')
########################################################################################
#Make composite plots of XXXXXXXXXXXXXXXXXXXXNNNCNNNNXXXXXXXXXXXXXXXXXXX
BWs <- c('../data/Tn5_XXXXXXXXXXXXXXXXXXXXNNNCNNNNXXXXXXXXXXXXXXXXXXX.bw')
#Reorder BWs to match Motifs:
compositelist <- vector('list', length(Motiflist))
for (i in 1:length(Motiflist)) {
  compositelist[[i]] = BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                        upstream = 100, downstream = 100,
                                        factor = names(Motiflist[i]),
                                        group = 'NNNCNNNN', ATAC = TRUE)
}
names(compositelist) = names(Motiflist)
Tn5_NNNCNNNN_compositelist = compositelist
save(Tn5_NNNCNNNN_compositelist, file = 'Tn5_NNNCNNNN_compositelist.Rdata')
########################################################################################
#Make composite plots of rule ensemble output
BWs <- c('Tn5_Rule_Ensemble_Scaled.bigWig')
#Reorder BWs to match Motifs:
compositelist <- vector('list', length(Motiflist))
for (i in 1:length(Motiflist)) {
  compositelist[[i]] = BED.query.bigWig(Motiflist[[i]], paste(BWs), paste(BWs), 
                                        upstream = 100, downstream = 100,
                                        factor = names(Motiflist[i]),
                                        group = 'Rule Ensemble', ATAC = TRUE)
}
names(compositelist) <- names(Motiflist)
Tn5_RE_compositelist = compositelist
save(Tn5_RE_compositelist, file = 'Tn5_RE_compositelist.Rdata')
#######################################################################################
#Plot difference b/t avg baseline and corrected...
#trim composites to 20bp around center
test_Tn5_RE = do.call(rbind, Tn5_RE_compositelist)
test_Tn5_RE$group = 'Rule Ensemble'
test_unscaled = do.call(rbind, Tn5_unscaled_compositelist)
test_unscaled$group = 'Unscaled'
test_NNNCNNNN = do.call(rbind, Tn5_NNNCNNNN_compositelist)
test_NNNCNNNN$group = 'seqOutBias'
###Subset out baseline values (motif +/-10) and single nucleotide (motif)
all_test_composites = rbind(test_Tn5_RE, test_unscaled, test_NNNCNNNN)

singlenuc_frame = NULL
singlenuc_store = NULL
for (i in 1:length(Motiflen)) {
  singlenuc_store = all_test_composites[which(all_test_composites$factor == 
                             names(Motiflen)[i] & all_test_composites$x < ((Motiflen[i]/2)+10) & 
                             all_test_composites$x > - ((Motiflen[i]/2)+10)),]
  singlenuc_frame = rbind(singlenuc_frame, singlenuc_store)
}
#Compute baseline and single nucleotide log2 ratios
singlenuc_log = NULL
for (i in 1:length(singlenuc_frame$factor)) {
  singlenuc_log[i] = log2(singlenuc_frame$est[i] / 0.006140407)
}
singlenuc_log = data.frame(singlenuc_log)
singlenuc_log[,2] = singlenuc_frame$factor
singlenuc_log[,3] = singlenuc_frame$group
colnames(singlenuc_log) = c('Difference', 'Factor', 'Treatment')
singlenuc_log$Treatment = factor(singlenuc_log$Treatment,
                                 levels = c('Unscaled', 'Rule Ensemble', 'seqOutBias'))
#Split main figure and supplemental
main_fig = c('AR', 'ESR1', 'MEF2A', 'MGA', 'NFATC3', 'NR2F2', 'REST', 'SP1', 'TEAD1')
supp_singlenuc_log = singlenuc_log[-c(which(singlenuc_log$Factor %in% main_fig )),]
supp_singlenuc_log$Treatment = factor(supp_singlenuc_log$Treatment,
                                      levels = c('Unscaled', 'seqOutBias', 'Rule Ensemble'))
singlenuc_log = singlenuc_log[which(singlenuc_log$Factor %in% main_fig ),]
singlenuc_log$Treatment = factor(singlenuc_log$Treatment,
                                 levels = c('Unscaled', 'seqOutBias', 'Rule Ensemble'))
#Plot each TF (singlenuc)
log2plot = bwplot(Difference ~ Treatment | Factor , data = singlenuc_log,
                  between=list(y=0.5, x = 0.5),
                  scales=list(x=list(draw=TRUE),rot = 45, alternating=c(1,1,1,1),cex=1.2,font=1),
                  #xlab = '',
                  ylim = c(-4.5, 3.5),
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
                    panel.violin.hack(..., col = c("#A1A3AB", "#FF0000", "#0000FF"),
                                      varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                    panel.stripplot(..., col='#54545380', do.out=FALSE,
                                    jitter.data=TRUE, amount = 0.2, pch = 16)
                    panel.bwplot(..., pch = '|', do.out = FALSE)
                    
                  })
pdf('Figure6C_Tn5_singlenuc_log2_comparison.pdf', useDingbats = FALSE, width=10.83, height=6)
trellis.par.set(box.umbrella = list(lty = 1, col="#93939300", lwd=2),
                box.rectangle = list(col = '#93939300', lwd=1.6),
                plot.symbol = list(col='#93939300', lwd=1.6, pch ='.'))
################################################################################
print(update(log2plot, layout=c(9,1)))
dev.off()
save(singlenuc_log, file = 'Figure6C_Tn5_singlenuc_log2_comparison.Rdata')
################################################################################
log2plot = bwplot(Difference ~ Treatment | Factor , data = supp_singlenuc_log,
                  between=list(y=0.5, x = 0.5),
                  scales=list(x=list(draw=TRUE),rot = 45, alternating=c(1,1,1,1),cex=1.2,font=1),
                  #xlab = '',
                  ylim = c(-4.5, 3.5),
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
                    panel.violin.hack(..., col = c("#A1A3AB", "#FF0000", "#0000FF"),
                                      varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                    panel.stripplot(..., col='#54545380',
                    do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                    panel.bwplot(..., pch = '|', do.out = FALSE)
                    
                  })

pdf('Supplemental_Figure6A_Tn5_singlenuc_log2_comparison.pdf',
    useDingbats = FALSE, width=10.83, height=6)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939300", lwd=2),
                box.rectangle = list(col = '#93939300', lwd=1.6),
                plot.symbol = list(col='#93939300', lwd=1.6, pch ='.'))

################################################################################
print(update(log2plot, layout=c(9,1)))
dev.off()
save(supp_singlenuc_log, file = 'Supplemental_Figure6C_Tn5_singlenuc_log2_comparison.Rdata')
##Combine rule ensemble, seqOutBias and unscaled compositelists for plotting
combined_compositelist <- vector(mode = 'list', length = 1)
for (i in 1:length(Tn5_unscaled_compositelist)) {
  combined_compositelist[[i]] <- rbind(Tn5_RE_compositelist[[i]],
                                       Tn5_unscaled_compositelist[[i]],
                                       Tn5_NNNCNNNN_compositelist[[i]])
}
names(combined_compositelist) <- names(Tn5_unscaled_compositelist)
#Separate out the TFs to be plotted (for 6C) and their motif lengths

figure6B_plot = do.call(rbind, combined_compositelist[c(4, 11, 16)])
figure6B_plot$group[which(figure6B_plot$group == 'NNNCNNNN')] = 'seqOutBias'
figure6B_plot$group = factor(figure6B_plot$group,
                                 levels = c('Rule Ensemble', 'seqOutBias', 'Unscaled'))
mlen <- Motiflen[c(4, 11, 16)]
mlen <- as.data.frame(mlen)
rownames(mlen) <- levels(as.factor(figure6B_plot$factor))


plot.composites(figure6B_plot, legend = TRUE, 
                pdf_name = 'Figure6B_Tn5_PRE_NNNCNNNN_comparison',
                figwidth = 8, figheight = 4,
                ylabel = '',
                xlabel = '',
                motifline = TRUE, Motiflen = mlen, layoutgrid = c(3,1),
                x_axis_range = -20:20, X_axis_ticks = seq(-20,20,10),
                Y_axis_ticks = seq(0,0.05,0.01), nYaxisdigits = 2,
                hline_val = 0.006140407, y_axis_range = seq(0,0.02,0.01),
                y_axis = FALSE, hline = TRUE, Y_ticks = FALSE, labsize = 0.85)

#Plot supplemental composites
supp_figure6B_plot = do.call(rbind, combined_compositelist[-c(4, 11, 16)])
supp_figure6B_plot$group[which(supp_figure6B_plot$group == 'NNNCNNNN')] = 'seqOutBias'
supp_figure6B_plot$group = factor(supp_figure6B_plot$group,
                             levels = c('Rule Ensemble', 'seqOutBias', 'Unscaled'))
mlen = Motiflen[-c(4, 11, 16)]
mlen = mlen[which(names(mlen) %in% test)]
mlen = as.data.frame(mlen)
rownames(mlen) = levels(as.factor(supp_figure6B_plot$factor))

plot.composites(supp_figure6B_plot, legend = TRUE, 
                pdf_name = 'Supplemental6B_DNasePN_PRE_NNNCNNNN_comparison',
                figwidth = 12, figheight = 10,
                ylabel = '',
                xlabel = '',
                motifline = TRUE, Motiflen = mlen, layoutgrid = c(5,3),
                x_axis_range = -20:20, X_axis_ticks = seq(-20,20,10),
                Y_ticks = FALSE,
                hline_val = 0.006140407, y_axis_range = seq(0,0.02,0.01),
                y_axis = FALSE, hline = TRUE, labsize = 0.85)

#Plot overlays:
################################################################################
################################################################################
#Plot seqOutBias overlay
seqOutBias_Tn5_RE_overlay = Tn5_NNNCNNNN_compositelist[which(names(Tn5_NNNCNNNN_compositelist) %in% test)]
seqOutBias_Tn5_RE_overlay = do.call(rbind, seqOutBias_Tn5_RE_overlay)
seqOutBias_Tn5_RE_overlay$correction = seqOutBias_Tn5_RE_overlay$factor
seqOutBias_Tn5_RE_overlay$factor = seqOutBias_Tn5_RE_overlay$group
seqOutBias_Tn5_RE_overlay$group = seqOutBias_Tn5_RE_overlay$correction
seqOutBias_Tn5_RE_overlay$factor = 'seqOutBias'
plot.composites(seqOutBias_Tn5_RE_overlay, legend = FALSE, 
                pdf_name = 'Figure6E_seqOutBias_Tn5_RE_overlay',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                hline_val = 0.006140407, Y_ticks = FALSE,col.lines = c("#8080804D"),
                hline = TRUE, hlinelwd = 10, hlinelty = 3,
                y_axis_range = seq(0,0.05,0.02), y_axis = TRUE)  
#Plot rule ensemble overlay
RE_Tn5_overlay = Tn5_RE_compositelist[which(names(Tn5_RE_compositelist) %in% test)]
RE_Tn5_overlay = do.call(rbind, RE_Tn5_overlay)
RE_Tn5_overlay$correction = RE_Tn5_overlay$factor
RE_Tn5_overlay$factor = RE_Tn5_overlay$group
RE_Tn5_overlay$group = RE_Tn5_overlay$correction
RE_Tn5_overlay$factor = 'Rule Ensemble'
plot.composites(RE_Tn5_overlay, legend = FALSE, 
                pdf_name = 'Figure6E_Tn5_RE_overlay',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                hline_val = 0.006140407, Y_ticks = FALSE,col.lines = c("#8080804D"),
                hline = TRUE, hlinelwd = 10, hlinelty = 3,
                y_axis_range = seq(0,0.05,0.02), y_axis = TRUE)  
#Plot unscaled overlay
unscaled_Tn5_overlay = Tn5_unscaled_compositelist[which(names(Tn5_unscaled_compositelist) %in% test)]
unscaled_Tn5_overlay = do.call(rbind, unscaled_Tn5_overlay)
unscaled_Tn5_overlay$correction = unscaled_Tn5_overlay$factor
unscaled_Tn5_overlay$factor = unscaled_Tn5_overlay$group
unscaled_Tn5_overlay$group = unscaled_Tn5_overlay$correction
unscaled_Tn5_overlay$factor = 'Unscaled'
plot.composites(unscaled_Tn5_overlay, legend = FALSE, 
                pdf_name = 'Figure6E_Tn5_unscaled_overlay',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                hline_val = 0.006140407, Y_ticks = FALSE,col.lines = c("#8080804D"),
                hline = TRUE, hlinelwd = 10, hlinelty = 3,
                y_axis_range = seq(0,0.05,0.02), y_axis = TRUE)  

#Plot Improvement of rule ensemble over seqOutBias
dot_frame = singlenuc_frame[which(singlenuc_frame$group == 'Unscaled')]
dot_frame = dot_frame[,-c(5)]
colnames(dot_frame) = c('Unscaled_x', 'Unscaled_est', 'Unscaled_factor', 'Unscaled_group')

dot_frame = cbind(dot_frame, singlenuc_frame[which(singlenuc_frame$group == 'Rule Ensemble'),1:4])
colnames(dot_frame)[5:8] = c('RE_x', 'RE_est', 'RE_factor', 'RE_group')
dot_frame = cbind(dot_frame, singlenuc_frame[which(singlenuc_frame$group == 'seqOutBias'),1:4])
colnames(dot_frame)[9:12] = c('sOB_x', 'sOB_est', 'sOB_factor', 'sOB_group')
identical(dot_frame$Unscaled_factor, dot_frame$RE_factor)

dot_frame = dot_frame[,-c(4,5,7,8,9,11,12)]
dot_frame$calc_avg = 0.006140407
dot_frame$RE_absdiff = abs(dot_frame$RE_est - dot_frame$calc_avg)
dot_frame$sOB_absdiff = abs(dot_frame$sOB_est - dot_frame$calc_avg)
rownames(dot_frame) = 1:nrow(dot_frame)

dot_frame$sOB_sub_RE = dot_frame$sOB_absdiff - dot_frame$RE_absdiff
dot_frame$plot_col  = 'RE Improvement over sOB'

length(which(dot_frame$sOB_absdiff > dot_frame$RE_absdiff))/ nrow(dot_frame)

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

