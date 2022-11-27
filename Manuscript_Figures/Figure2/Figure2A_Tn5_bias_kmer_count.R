source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
library(lattice)
options(scipen=10000)

setwd('../output')

load('Figure2A_Upstream_OE_kmers.Rdata')

load('220801_Tn5_kmer_count_correction.Rdata')

mertable$upOE = mertable$upOE*Tn5_kmer_count_correction
mertable$downOE = mertable$downOE*Tn5_kmer_count_correction

pdf(file = 'Figure2A_Tn5_upstream_OE_kmers.pdf', width=20, height=5)
barchart(upOE~kmer, mertable,
         ylab = list(expression(bold(paste(frac(Observed, Expected), " 3-mers" ))), cex = 3),
         scales=list(x=list(cex=2.2, rot = 90, fontfamily = "mono"), y=list(cex=2.5)),
         ylim = c(0,1.8),
         col = 'black',
         par.settings = list(axis.line = list(col = 0)),
         panel=function(...) {
           panel.barchart(...)
           panel.abline(h=1, col = 'red', lwd = 4, lty = 2)
         })
dev.off()
