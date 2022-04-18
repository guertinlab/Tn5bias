library(lattice)
library(ggplot2)
################################################################################
################################################################################
################################################################################
################################################################################
plot.composites = function(dat, ylabel = '', pdf_name = 'PLEASE_SET_FILE_NAME',
         xlabel = '', striplabel = TRUE, legend = TRUE,
         motifline = FALSE, Motiflen = 10,
         figwidth = 2.5, figheight=3,
         indexlist = NULL, layoutgrid = NULL,
         col.lines = c("#0000FF", "#FF0000", "#00000090", 
                       rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                       rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), 
                       rgb(1/2,1/2,0,1/2)), 
         fill.poly = c(rgb(0,0,1,1/4), 
                       rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4), 
                       rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
  require(lattice)
  pdf(paste(pdf_name, '.pdf', sep = ''), width= figwidth, height= figheight) 
  print(xyplot(est ~ x|factor, group = group, data = dat, strip = striplabel,
               type = 'l', as.table = TRUE,
               scales=list(x=list(cex=1.5,relation = "free", axs ="i"), 
                           y =list(cex=1.5, relation="free", tick.number=4)),
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
               par.strip.text=list(cex=2, font=1, col='black'),
               ylim = c(0, 0.025),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               lwd=2,
               ylab = list(label = paste(ylabel), cex =1.75),
               xlab = list(label = paste(xlabel), cex =1.75),
               index.cond = indexlist,
               layout = layoutgrid,
               panel = function(x, y, ...) {
                 panel.xyplot(x, y, ...)
                 #panel.abline(h = 0, lty =1, lwd = 1.0, col = '#A9A9A932')
                 if (motifline == TRUE) 
                 {panel.abline(v = Motiflen/2, lty = 2, col = "red")} else{}
                 if (motifline == TRUE) 
                 {panel.abline(v = -Motiflen/2, lty = 2, col = "red")} else{}
               }
  ))
  dev.off()
}
################################################################################
################################################################################
train = c('AR', 'ASCL1', 'ATF3', 'CEBPB', 'CLOCK', 'CTCF', 'DLX2', 'DUX4', 'EGR1', 'ELF1', 'ESR1', 'FOXA1', 'FOXK2', 'GLIS1', 'HSF1', 'JUN', 'LEF1', 'MEIS2', 'MLX', 'MYC', 'NR2F2', 'SPIB', 'SRY', 'STAT1', 'TEAD1', 'TGIF1')
test = c('E2F1', 'EBF1', 'FERD3L', 'FOS', 'GATA2', 'HOXC12', 'IRF1', 'MAX', 'MEF2A', 'MGA', 'NFATC3', 'POU3F1', 'PPARG', 'REST', 'RUNX1', 'Six3', 'SP1', 'USF1')



################################################################################
load('Figure5A_Tn5_unscaled_composites_overlay.Rdata')

plot.composites(Tn5_unscaled_composites, legend = FALSE, 
                  pdf_name = 'Figure5A_Tn5_unscaled_overlay',
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8)  

  
################################################################################
load('Figure5A_Tn5_NNNCNN_composites_overlay.Rdata')
plot.composites(Tn5_NNNCNN_composites, legend = FALSE, 
                  pdf_name = 'Figure5A_Tn5_NNNCNN_overlay',
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8)  


################################################################################
load('Figure5A_Tn5_PRE_composites_overlay.Rdata')  
plot.composites(Tn5_PRE_composites, legend = FALSE, 
                  pdf_name = 'Figure5A_Tn5_PRE_overlay',
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8)  


################################################################################
load('Figure5B_Baseline_quantification.Rdata')
Maskvarplot <- function(){
    ggplot(plotvar, aes(x=factor(Method), y=Avg)) +
      xlab('') + ylab('Tn5 Cut Frequency') +
      geom_boxplot(fill= '#6495ED4D', outlier.size = 0.8, lwd=1.25, fatten = 1) +
      geom_point(size = 1, shape = 20) +
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 0.650, hjust=0.75, size = 16, face = 'bold', color = 'black'), 
            axis.title.y=element_text(size=16, color = 'black', face = 'bold'),
            axis.title.x = element_text(vjust = -1),
            axis.text.y = element_text(size = 16, color = 'black', face = 'bold'),
            axis.line = element_line(colour = 'black', size = 1)) 
  }
  
  pdf("Figure5B_Tn5_unscaled_PRE_NNNCNN_basemean.pdf", width=6, height=7)
  Maskvarplot()
  dev.off()
  
  