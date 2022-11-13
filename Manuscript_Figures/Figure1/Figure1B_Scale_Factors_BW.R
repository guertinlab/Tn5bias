library(lattice)
setwd('../output')
#Load in  scale factor files
load('../data/Figure1B_nuclease_scalefactors_list.Rdata')

for (i in 1:length(Figure1C_list)) {
  Figure1C_list[[i]][,14] = as.factor(Figure1C_list[[i]][,14])
}
############################################################################################################

bw_theme <- trellis.par.get()
bw_theme$box.dot$pch <- "|"
bw_theme$box.rectangle$fill <- "light blue"
bw_theme$box.umbrella$lty <- 1
bw_theme$box.umbrella$lwd <- 3
bw_theme$box.umbrella$col <- "black"
bw_theme$box.rectangle$col <- "black"
bw_theme$box.rectangle$lwd <- 3


for (i in 1:length(Figure1C_list)) {
pdf(paste(names(Figure1C_list)[i], '_maskpositions.pdf', sep = ''), width=26, height=7)
print(bwplot(log2(1/plusscalefact)~xaxis, Figure1C_list[[i]],
       scales=list(x=list(rot=90), y=list(rot=0), relation="free", cex=3.5, font=1),
       ylab="",
       main="",
       do.out=FALSE,
       par.settings = bw_theme,
       panel=function(...){
         panel.bwplot(...)
         panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, cex = 0.33, pch = 16, alpha = 0.5)
       }))
dev.off()
}
