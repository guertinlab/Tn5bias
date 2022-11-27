library(lattice)

setwd('../output')
load('../data/Figure2B_CAG_correction.Rdata')

##Plot.composites takes the composite.lattice object and creates a plot for 
#this data, while also allowing a mapping of the original motif width onto the
#plot by specifying a TFlen value (the default is 10). 
plot.composites <- function(dat, ylabel = '', pdf_name = 'PLEASE_SET_FILE_NAME',
                            xlabel = '', striplabel = TRUE, legend = TRUE,
                            motifline = FALSE, Motiflen = 10,
                            figwidth = 2.5, figheight=3,
                            indexlist = NULL, layoutgrid = NULL,
                            col.lines = c("#0000FF", "#0000FF", "#0000FF", "#a1a3ab"), 
                            fill.poly = c(rgb(0,0,1,1/4), 
                                          rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4), 
                                          rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
  require(lattice)
  pdf(paste(pdf_name, '.pdf', sep = ''), width= figwidth, height= figheight) 
  print(xyplot(est ~ x|factor, group = group, data = dat, strip = striplabel,
               type = 'l', as.table = TRUE,
               scales=list(x=list(at=seq(-15,15,5), cex=1.4, relation = "free", axs ="i", rot = 90), 
                           y =list(cex=1.4, relation="free", tick.number=2)),
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
               cex.axis=1.4,
               par.strip.text=list(cex=1.2, font=1, col='black'),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               lwd=3,
               ylab = list(label = paste(ylabel), cex =1.4),
               xlab = list(label = paste(xlabel), cex =1.4),
               xlim = c(-15,15),
               ylim = c(0,0.0175),
               index.cond = indexlist,
               layout = layoutgrid,
               panel = function(x, y, ...) {
                 panel.xyplot(x, y, ...)
                 #panel.abline(h = 0, lty =1, lwd = 1.0, col = '#A9A9A932')
                 level = dimnames(trellis.last.object())[["factor"]][packet.number()]
                 if (motifline == TRUE) 
                 {panel.abline(v = Motiflen[rownames(Motiflen)==level,]/2, lty = 2, col = "red")} else{}
                 if (motifline == TRUE) 
                 {panel.abline(v = -Motiflen[rownames(Motiflen)==level,]/2, lty = 2, col = "red")} else{}
               }
  ))
  dev.off()
}






plot.composites(Tn5_CAG_peak_direction_composite, legend = FALSE, 
                pdf_name = 'Figure2B_CAG_direction_maskcompare',
                ylabel = 'Cut Frequency',
                xlabel = 'Distance from CAG Center',
                indexlist = list(c(3,2,1)),
                layoutgrid = c(3,1),
                figwidth = 8, figheight=5,
                motifline = FALSE)
