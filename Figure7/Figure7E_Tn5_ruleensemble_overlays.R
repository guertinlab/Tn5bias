library(lattice)
library(data.table)
library(grid)
setwd('../output')
################################################################################

load('../data/Figure7E_Tn5_unscaled_overlay.Rdata')
load('../data/Figure7E_Tn5_RuleEnsemble_overlay.Rdata')
load('../data/Figure7E_Tn5_seqOutBias_overlay.Rdata')
################################################################################
################################################################################
plot.composites = function(dat, ylabel = '', x_axis_range = min(dat$x):max(dat$x), pdf_name = 'PLEASE_SET_FILE_NAME',
                           xlabel = '', striplabel = TRUE, legend = TRUE,
                           motifline = FALSE, Motiflen = 10,
                           figwidth = 2.5, figheight=3,
                           indexlist = NULL, layoutgrid = NULL,
                           col.lines = c('#8080804D'), 
                           fill.poly = c(rgb(0,0,1,1/4), 
                                         rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4), 
                                         rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
  require(lattice)
  dat = dat[which(dat$x %in% x_axis_range),]
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
               par.strip.text=list(cex=2.5, font=1, col='black'),
               ylim = c(0, 0.05),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               lwd=2,
               ylab = list(label = paste(ylabel), cex =3.00),
               xlab = list(label = paste(xlabel), cex =3.00),
               index.cond = indexlist,
               layout = layoutgrid,
               panel = function(x, y, ...) {
                 panel.xyplot(x, y, ...)
                 panel.abline(h = 0.006140407, lty =3, lwd = 10.0, col = 'red')
                 
               }
  ))
  dev.off()
}


################################################################################
################################################################################

plot.composites(OLD_Tn5_RE_overlay, legend = FALSE, 
                pdf_name = 'Figure7E_seqOutBias_NNNCNN_Tn5_RE_overlay',
                ylabel = 'Cut Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8)  


plot.composites(RE_Tn5_overlay, legend = FALSE, 
                pdf_name = 'Figure7E_Tn5_RE_overlay',
                ylabel = 'Cut Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8)  


plot.composites(unscaled_Tn5_overlay, legend = FALSE, 
                pdf_name = 'Figure7E_unscaled_Tn5_overlay',
                ylabel = 'Cut Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8)  



