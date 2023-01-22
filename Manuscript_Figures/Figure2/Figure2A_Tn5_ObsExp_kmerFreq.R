source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
library(data.table)
library(lattice)
library(stringr)
options(scipen=10000)
#
load('Figure2A_Tn5_upstream_OE_kmers.Rdata')
#Plot corrected upstream observed/expected k-mer counts
pdf(file = 'Figure2A_Tn5_upstream_OE_kmers.pdf', width=20, height=5)
barchart(CorrupOE~kmer, mertable,
         ylab = list(expression("log"[2]~frac(Observed, Expected)~" 3-mers" ), cex = 3),
         scales=list(x=list(cex=2.2, rot = 90, fontfamily = "mono"), y=list(cex=2.5)),
         ylim = c(-1,1),
         col = 'black',
         par.settings = list(axis.line = list(col = 0)),
         panel=function(...) {
           panel.barchart(...)
           panel.abline(h=0, col = 'red', lwd = 4, lty = 2)
         })
dev.off()
