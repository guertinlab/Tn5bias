library(ggplot2)
library(lattice)
setwd('../output')
################################################################################
################################################################################
################################################################################
################################################################################
plot.composites = function(dat, ylabel = '', pdf_name = 'PLEASE_SET_FILE_NAME',
                           xlabel = '', striplabel = TRUE, legend = TRUE,
                           motifline = FALSE, Motiflen = 10,
                           figwidth = 2.5, figheight=3,
                           indexlist = NULL, layoutgrid = NULL,
                           col.lines = c('#8080804D'), 
                           fill.poly = c(rgb(0,0,1,1/4), 
                                         rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4), 
                                         rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
  require(lattice)
  pdf(paste(pdf_name, '.pdf', sep = ''), width= figwidth, height= figheight) 
  print(xyplot(est ~ x|factor, group = group, data = dat, strip = striplabel,
               type = 'l', as.table = TRUE,
               scales=list(x=list(cex=2.5,relation = "free", axs ="i"), 
                           y =list(cex=2.5, relation="free", tick.number=4)),
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
               ylim = c(0, 0.05),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               lwd=2,
               ylab = list(label = paste(ylabel), cex =3.0),
               xlab = list(label = paste(xlabel), cex =3.0),
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
load('../data/Figure7A_DNase_unscaled_composites_overlay.Rdata')

plot.composites(DNase_unscaled_composites_overlay, legend = FALSE, 
                  pdf_name = 'Figure7A_DNase_unscaled_overlay_mk4',
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8)


################################################################################
load('../data/Figure7A_Tn5_unscaled_composites_overlay.Rdata')  

plot.composites(Tn5_unscaled_composites_overlay, legend = FALSE, 
                  pdf_name = 'Figure7A_Tn5_unscaled_overlay_mk4',
                  ylabel = 'Cut Frequency',
                  xlabel = 'Distance from Motif Center',
                  motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8)  

  

  
