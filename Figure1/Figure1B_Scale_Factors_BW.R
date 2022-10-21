library(ggplot2)
library(reshape2)
library(data.table)

#Load in  scale factor files
load('Figure1B_nuclease_scalefactors_list.Rdata')

#Make lattice plot data.tables:
Benzonase_plot = Figure1C_list[[1]]
Cyanase_plot = Figure1C_list[[2]]
MNase_plot = Figure1C_list[[3]]
DNase_plot = Figure1C_list[[4]]
Tn5_plot = Figure1C_list[[5]]


for (i in 1:(length(Figure1C_list)-1)) {
  pdf(paste(names(Figure1C_list)[i], '_maskpositions.pdf', sep = ''), width=23, height=7)
  par(mar=c(10,4,4,4))
  boxplot(log2((1/plusscalefact))~xaxis,data=Figure1C_list[[i]], las=1, pch = 16, outcex=0.33, pars  =  list(xaxt = "n"),
          ylab = '',
          xlab = '',
          whisklty = 1,
          whisklwd = 4,
          staplelwd = 4,
          ylim = c(-6,4),
          boxcol = 'light blue',
          col = 'light blue',
          frame.plot = FALSE,
          par(cex.axis=3))
  axis(1, at=seq(1, 32, by=1), labels = FALSE, tck = -0.05, lwd.ticks = 2)
  axis(2, lwd.ticks = 2, labels = FALSE)
  text(x = seq(1, 32, by=1)+0.25, par("usr")[3] - 0.75, labels = seq(-15.5, 15.5, by=1), srt = 90, pos = 2, xpd = TRUE, cex = 3)
  box(bty="l")
  dev.off()
}




pdf('Tn5_scalefactors_maskpositions.pdf', width=23, height=7)
par(mar=c(10,4,4,4))
boxplot(log2((1/plusscalefact))~xaxis,data=Figure1C_list[[5]], las=1, pch = 16, outcex=0.33, pars  =  list(xaxt = "n"),
        ylab = '',
        xlab = '',
        whisklty = 1,
        whisklwd = 4,
        staplelwd = 4,
        ylim = c(-6,4),
        boxcol = 'light blue',
        col = 'light blue',
        frame.plot = FALSE,
        par(cex.axis=3))
axis(1, at=seq(1, 31, by=1), labels = FALSE, tck = -0.05, lwd.ticks = 2)
axis(2, lwd.ticks = 2, labels = FALSE)
text(x = seq(1, 31, by=1)+0.25, par("usr")[3] - 0.75, labels = seq(-15, 15, by=1), srt = 90, pos = 2, xpd = TRUE, cex = 3)
box(bty="l")
dev.off()


