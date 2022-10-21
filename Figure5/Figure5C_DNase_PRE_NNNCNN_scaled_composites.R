library(ggplot2)
setwd('../output')
################################################################################
################################################################################
################################################################################
################################################################################
##Plot.composites takes the composite.lattice object and creates a plot for 
#this data, while also allowing a mapping of the original motif width onto the
#plot by specifying a TFlen value (the default is 10). 
plot.composites <- function(dat, ylabel = '', pdf_name = 'PLEASE_SET_FILE_NAME',
                            xlabel = '', striplabel = TRUE, legend = TRUE,
                            motifline = FALSE, Motiflen = 10,
                            figwidth = 2.5, figheight=3,
                            indexlist = NULL, layoutgrid = NULL,
                            col.lines = c("#0000FF", "#FF0000",  
                                          "#00ffaf", rgb(0,0,0,1/2),  
                                          rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), 
                                          rgb(1/2,1/2,0,1/2)), 
                            fill.poly = c(rgb(0,0,1,1/4), 
                                          rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4), 
                                          rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
  require(lattice)
  pdf(paste(pdf_name, '.pdf', sep = ''), width= figwidth, height= figheight) 
  print(xyplot(est ~ x|factor, group = group, data = dat, strip = striplabel,
               type = 'l', as.table = TRUE,
               scales=list(x=list(cex=0.85,relation = "free", axs ="i"), 
                           y =list(cex=0.85, relation="free", tick.number=4)),
               col = col.lines,
               auto.key = if (legend == TRUE) 
               {list(points=F, lines=T, cex=0.8)} else{},
               par.settings = list(strip.background=list(col="#00000000"),
                                   strip.border = list(col = 'transparent'),
                                   superpose.symbol = list(pch = c(16),
                                                           col=col.lines, 
                                                           cex =0.5), 
                                   superpose.line = list(col = col.lines, 
                                                         lwd=c(2), 
                                                         lty = c(1))),
               cex.axis=1.0,
               par.strip.text=list(cex=1.2, font=1, col='black'),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               lwd=2,
               ylab = list(label = paste(ylabel), cex =1.2),
               xlab = list(label = paste(xlabel), cex =1.2),
               index.cond = indexlist,
               layout = layoutgrid,
               panel = function(x, y, ...) {
                 panel.xyplot(x, y, ...)
                 #panel.abline(h = 0, lty =1, lwd = 1.0, col = '#A9A9A932')
                 level = dimnames(trellis.last.object())[["factor"]][packet.number()]
                 if (motifline == TRUE) 
                 {panel.abline(v = ceiling(Motiflen[rownames(Motiflen)==level,]/2), lty = 2, col = "red")} else{}
                 if (motifline == TRUE) 
                 {panel.abline(v = ceiling(-Motiflen[rownames(Motiflen)==level,]/2), lty = 2, col = "red")} else{}
               }
  ))
  dev.off()
}


load('../data/Figure5C_DNase_PRE_NNNCNN_unscaled_compositelist_mk2.Rdata')
load('../data/Figure5C_Motif_lengths.Rdata')

test = c('AR', 'CEBPB', 'E2F1', 'ESR1', 'FERD3L', 'HOXC12', 'HSF1', 'MAX', 'MEF2A', 'MGA', 'NFATC3', 'NR2F2', 'POU3F1', 'REST', 'Six3', 'SP1', 'USF1', 'TEAD1')
test_Motiflen = Motiflen[which(names(Motiflen) %in% test)]
Motiflen = test_Motiflen
#trim composites to 20bp around center
for (i in 1:length(combined_compositelist)) {
  combined_compositelist[[i]] = combined_compositelist[[i]][which(combined_compositelist[[i]]$x >= -20.5 & combined_compositelist[[i]]$x <= 20.5),]
}

#Separate out the TFs to be plotted and their motif lengths
Figure5C_plot = do.call(rbind, combined_compositelist[c(2, 7, 18)])
Figure5C_plot[which(Figure5C_plot$group == 'Rules Ensemble'),4] = 'Rule Ensemble'
Figure5C_plot[which(Figure5C_plot$group == 'NNNCNN'),4] = 'seqOutBias'
Figure5C_plot[which(Figure5C_plot$group == 'unscaled'),4] = 'Unscaled'
mlen <- Motiflen[c(2, 7, 18)]
mlen <- as.data.frame(mlen)
rownames(mlen) <- levels(as.factor(Figure5C_plot$factor))



plot.composites(Figure5C_plot, legend = TRUE, 
                  pdf_name = 'Figure5C_DNasePN_PRE_NNNCNN_comparison_mk1',
                  figwidth = 8, figheight = 4,
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = TRUE, Motiflen = mlen, layoutgrid = c(3,1))

