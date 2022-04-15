library(ggplot2)
library(ggseqlogo)

plot.seqlogo.func <- function(x, outfile = "ATAC-kmer_optimization_all_test.pdf") {
  w =  0.663 + (ncol(x) + 1)*0.018 + (ncol(x)+2)* .336
  logoseq = ggseqlogo(x,  facet = "wrap", font = 'helvetica_bold') 
  pdf(outfile, useDingbats=FALSE, width=w, height=2.695)
  print(logoseq + scale_x_continuous(breaks =  c(seq(1,ncol(x),1)), 
                               labels = c(seq(-(ncol(x) %/% 2),(ncol(x) %/% 2),by=1))))
  dev.off()
}

###################################################################################
setwd('../git')




###################################################################################
load('../git/Benzonase_sepcat_bias_pswm_31win.Rdata')
plot.seqlogo.func(pswm.Benzonase_sepcat.trans, outfile='pswm_Benzonase_sepcat.pdf')


###################################################################################
load('../git/Cyanase_sepcat_bias_pswm_31win.Rdata')
plot.seqlogo.func(pswm.Cyanase_sepcat.trans, outfile='pswm_Cyanase_sepcat.pdf')



###################################################################################
load('../git/DNase_sepcat_bias_pswm_31win.Rdata')
plot.seqlogo.func(pswm.DNase_sepcat.trans, outfile='pswm_DNase_sepcat.pdf')



###################################################################################
load('../git/MNase_sepcat_bias_pswm_31win.Rdata')
plot.seqlogo.func(pswm.MNase_sepcat.trans, outfile='pswm_MNase_sepcat.pdf')




###################################################################################
plot.seqlogo.func <- function(x, outfile = "ATAC-kmer_optimization_all_test.pdf") {
  w =  0.663 + (ncol(x) + 1)*0.018 + (ncol(x)+2)* .336
  logoseq = ggseqlogo(x,  facet = "wrap", font = 'helvetica_bold') 
  pdf(outfile, useDingbats=FALSE, width=w, height=2.695)
  print(logoseq + scale_x_continuous(breaks =  c(seq(1,ncol(x),1)), 
                                     labels = c(seq(-(ncol(x) %/% 2),(ncol(x) %/% 2),by=1))) + 
                  xlab("Distance from Centrally Recognized Base") +
                  theme(axis.title.x = element_text(size=14, face="bold", vjust = -1)))
  dev.off()
}


load('../git/tn5_sepcat_bias_pswm.Rdata')
plot.seqlogo.func(pswm.Tn5_sepcat.trans, outfile='pswm_Tn5_sepcat.pdf')



