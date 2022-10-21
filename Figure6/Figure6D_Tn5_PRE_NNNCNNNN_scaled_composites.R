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
                                          "#00ffaf", rgb(0,0,0,1/2),"#3d7a72", "#96447f","#ff908a", rgb(0,0,0,1/2),  
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
                 panel.abline(h = 0.006140407, lty = 2, col = "red", lwd = 2)
               }
  ))
  dev.off()
}


load('../data/Figure6E_Tn5_RE_NNNCNNNN_unscaled_compositelist_mk1.Rdata')
load('../data/Figure6E_Motif_lengths.Rdata')

#trim composites to 20bp around center
for (i in 1:length(combined_compositelist)) {
  combined_compositelist[[i]] = combined_compositelist[[i]][which(combined_compositelist[[i]]$x >= -20 & combined_compositelist[[i]]$x <= 20),]
}

#Separate out the TFs to be plotted and their motif lengths
figure6E_plot = do.call(rbind, combined_compositelist[c(13, 20, 38)])
figure6E_plot[which(figure6E_plot$group == 'NNNCNNNN'),4] = 'seqOutBias'
figure6E_plot[which(figure6E_plot$group == 'unscaled'),4] = 'Unscaled'
mlen <- Motiflen[c(13, 20, 38)]
mlen <- as.data.frame(mlen)
rownames(mlen) <- levels(as.factor(figure6E_plot$factor))



plot.composites(figure6E_plot, legend = TRUE, 
                  pdf_name = 'Figure6C_Tn5_PRE_NNNCNN_comparison_mk1',
                  figwidth = 8, figheight = 4,
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = TRUE, Motiflen = mlen, layoutgrid = c(3,1))

