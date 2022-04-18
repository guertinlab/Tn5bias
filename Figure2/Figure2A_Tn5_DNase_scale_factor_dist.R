library(ggplot2)




load('Tn5_5mer_scalefactors.Rdata')

ATAC.plot.new <- function(){
  ggplot(Tn5_5mer_scalefactors, aes(x=factor(xaxis), y=value)) +
    xlab('Scale Factor Position Relative to Central Base') + ylab(expression(log[2]*(frac(1,'scale factor')))) + ylim(c(-6, 4)) +
    geom_boxplot(fill= '#6495ED', outlier.size = 0.8) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.0, hjust=0.75, size = 9, face = 'bold', color = 'black'), 
          axis.title.y=element_text(size=16, face = 'bold'),
          axis.title.x = element_text(vjust = -1.75),
          axis.text.y = element_text(angle = 0, vjust = 0.0, hjust=0.0, size = 12, face = 'bold', color = 'black')) 
}

pdf("ATAC_Allscalingvalues.pdf", width=9, height=7)
ATAC.plot.new()
dev.off()




######################################################################################################################################


load('DNase_5mer_scalefactors.Rdata')

DNase.plot.new <- function(){
  ggplot(DNase_5mer_scalefactors, aes(x=factor(xaxis), y=value)) +
    xlab('Scale Factor Position Relative to Cutsite') + ylab(expression(log[2]*(frac(1,'scale factor')))) + ylim(c(-6, 4)) +
    geom_boxplot(fill= '#6495ED', outlier.size = 0.8) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.0, hjust=0.75, size = 9, face = 'bold', color = 'black'), 
          axis.title.y=element_text(size=16, face = 'bold'),
          axis.title.x = element_text(vjust = -1.75),
          axis.text.y = element_text(angle = 0, vjust = 0.0, hjust=0.0, size = 12, face = 'bold', color = 'black'))
}

pdf("DNase_Allscalingvalues.pdf", width=9, height=7)
DNase.plot.new()
dev.off()
