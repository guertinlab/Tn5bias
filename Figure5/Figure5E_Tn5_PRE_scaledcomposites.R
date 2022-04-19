library(lattice)
source('composite_functions.R')
################################################################################
################################################################################
################################################################################
################################################################################

#
train = c('AR', 'ASCL1', 'ATF3', 'CEBPB', 'CLOCK', 'CTCF', 'DLX2', 'DUX4', 'EGR1', 'ELF1', 'ESR1', 'FOXA1', 'FOXK2', 'GLIS1', 'HSF1', 'JUN', 'LEF1', 'MEIS2', 'MLX', 'MYC', 'NR2F2', 'SPIB', 'SRY', 'STAT1', 'TEAD1', 'TGIF1')
test = c('E2F1', 'EBF1', 'FERD3L', 'FOS', 'GATA2', 'HOXC12', 'IRF1', 'MAX', 'MEF2A', 'MGA', 'NFATC3', 'POU3F1', 'PPARG', 'REST', 'RUNX1', 'Six3', 'SP1', 'USF1')



load('Figure5E_motiflengths.Rdata')

load('Figure5E_Tn5_PRE_NNNCNN_unscaled_composites.Rdata')
                                        #trim composites to 20bp around center
x = combined_compositelist[[1]]
for (i in 2:length(combined_compositelist)) {
  x = rbind(x, combined_compositelist[[i]])
}

data.composite = x

save(data.composite, file = "data.composite.Rdata")


plot.composites <- function(dat, ylabel = '', pdf_name = 'PLEASE_SET_FILE_NAME',
                            xlabel = '', striplabel = TRUE, legend = TRUE,
                            motifline = FALSE, Motiflen = 10,
                            figwidth = 2.5, figheight=3,
                            indexlist = NULL, layoutgrid = NULL, xlim = c(-21, 21),
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
                 type = 'l', as.table = TRUE, xlim = xlim,
                 scales=list(x=list(cex=0.8,relation = "free", axs ="i"), 
                             y =list(cex=0.8, relation="free", tick.number=4)),
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
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 lwd=2,
                 ylab = list(label = paste(ylabel), cex =0.8),
                 xlab = list(label = paste(xlabel), cex =0.8),
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


x = data.composite
x[,4] = gsub("NNNCNN", "seqOutBias", x[,4])
x[,4] = gsub("PRE", "Rule Ensemble", x[,4])

x[,4] <- as.factor(x[,4])


for (i in c('Rule Ensemble', 'seqOutBias', 'Unscaled' )) {
  x[,4] = relevel(x[,4],ref=i)
}


plot.composites(x,
                  pdf_name = 'composite',
                  ylabel = 'Relative Cleavage', indexlist = list(c(10,13,35,44)), 
                  xlim = c(-21, 21),  figwidth = 11, figheight=8,
                  xlabel = 'Distance from Motif Center',
  )


