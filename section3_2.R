library(parallel)
library(devtools)
library(lattice)
devtools::install_github("andrelmartins/bigWig", subdir="bigWig")
source("https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_hcsearch.R")

load.sites <- function(filenames) {
  lapply(filenames, function(filename) {
    fimo = read.table(filename, skip = 1, sep ='\t')
    len.vec.values = length(fimo[,7])
    if (length(fimo[,7]) > 3000000 ){
      subset.len = len.vec.values - 300000
      sort.sub = sort(fimo[,7], partial=subset.len)[subset.len]
      top.index = which(fimo[,7] > sort.sub)
      fimo.small = fimo[top.index,]
      bed6 = fimo.small[, c(3,4,5,2,7,6)]
    } else {
      bed6 = fimo[, c(3,4,5,2,7,6)]
    }
    bed6
  }) 
}

setwd('rep1_HC')

#these are the top 13 positions identified by hill climbing:
hc.atac.cutmasks.rep1 = c('XXXXXXXXXXXXXXXXXXXXX',
                          'XXXXXXNXXXXXXXXXXXXXX',
                          'XXXXXXNXXXXXXXNXXXXXX',
                          'XXXXXXNXXXXNXXNXXXXXX',
                          'XXXXXXNXXNXNXXNXXXXXX',
                          'XXXXNXNXXNXNXXNXXXXXX',
                          'XNXXNXNXXNXNXXNXXXXXX',
                          'XNXXNXNXXNXNXXNXNXXXX',
                          'XNXXNXNNXNXNXXNXNXXXX',
                          'XNXXNXNNXNXNNXNXNXXXX',
                          'XNXXNXNNXNXNNXNXNXXNX',
                          'XNXXNXNNXNNNNXNXNXXNX',
                          'XNXXNXNNXNNNNNNXNXXNX',
                          'XNNXNXNNXNNNNNNXNXXNX')
		

sites = load.sites(c("../rep1_ATAC_bias_pe1_plus_00005_fimo.txt", "../rep1_ATAC_bias_pe1_minus_00005_fimo.txt",
                     "../rep1_ATAC_bias_pe2_plus_00005_fimo.txt", "../rep1_ATAC_bias_pe2_minus_00005_fimo.txt"))

hc.scores.atac.rep1 = mclapply(hc.atac.cutmasks.rep1, function(cutmask) {
	bw.plus = load.bigWig(paste0('runhc_21merrep1_', cutmask, '.bigWig'))
	bw.minus = load.bigWig(paste0('runhc_21merrep1_', cutmask, '.bigWig'))
	eval.cutmask(sites, bw.plus, bw.minus)
}, mc.cores = 4)


save(hc.scores.atac.rep1,
	file = 'hc_scores_atac_kmer_optimization_21mer_rep1.Rdata')

hc.atac.cutmasks.rep1 = factor(hc.atac.cutmasks.rep1, levels=hc.atac.cutmasks.rep1)

pdf("HC_scores_atac_kmer_optimization_21mer_rep1.pdf", useDingbats=FALSE, width=6, height=6)
dotplot(as.numeric(hc.scores.atac.rep1) ~ hc.atac.cutmasks.rep1,
	pch = 19,
	cex =1,
	col = 'black',
	main = "Hill Climbing derived k-mer masks - rep1",
	xlab = 'Masked Positions',
	ylim = c(0, 150000),
	scales=list(x=list(rot=65)),
	ylab = expression(paste(Sigma, ' SDs between PSWM positions')))
#for each set of TF PSWMs we sum the intensity of signal at each position, 
#then we take the standard deviation between positions.
#The final metric is a sum of these standard deviations.
dev.off()
