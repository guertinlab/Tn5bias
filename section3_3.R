source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
library(lattice)

parse.fimo <- function(file) {
    fimo.data = read.table(file, skip = 1, sep ='\t')
    res = fimo.data[, c(3,4,5,8,9,6, 2)]
    colnames(res) = c('chr', 'start', 'end', 'score', 'pval', 'strand', 'motif')
    return(res)
}


setwd('Directory Containing .bigWig files')

system('mkdir final')
setwd('final')
system('cp ../runhc_21merrep1_XXXXXXXXXXXXXXXXXXXXX.bigWig ./')
system('cp ../runhc_21merrep1_XXXXXXNXXNXNXXNXXXXXX.bigWig ./')
system('cp ../runhc_21merrep1_XNXXNXNXXNXNXXNXNXXXX.bigWig ./')
system('cp ../runhc_21merrep1_XNXXNXNNXNNNNXNXNXXNX.bigWig ./')


all.composites.ATAC.bias.rep1.21mer = cycle.fimo.new.not.hotspots(path.dir.fimo = 'Path to FIMO .txt files', 
                   path.dir.bigWig = './final', 
                   window = 30, exp = 'ATAC')

save(all.composites.ATAC.bias.rep1.21mer, file = "all.composites.ATAC.bias.rep1.21mer.Rdata")
