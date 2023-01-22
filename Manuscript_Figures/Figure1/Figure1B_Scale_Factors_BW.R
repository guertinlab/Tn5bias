library(lattice)
#Load in  scale factor files
load('Figure1B_scalefactor_list.Rdata')

for (i in 1:length(Figure1B_list)) {
  Figure1B_list[[i]][,2] = as.factor(Figure1B_list[[i]][,2])
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

for (i in 1:length(Figure1B_list)) {
  pdf(paste(names(Figure1B_list)[i], '_maskpositions.pdf', sep = ''), width=26, height=7)
  print(bwplot(log2(1/minusscalefact)~xaxis, Figure1B_list[[i]],
               scales=list(x=list(rot=90), y=list(rot=0), relation="free", cex=3.5, font=1),
               ylab="",
               main="",
               do.out=FALSE,
               par.settings = bw_theme,
               panel=function(...){
                 panel.bwplot(...)
               }))
  dev.off()
}
