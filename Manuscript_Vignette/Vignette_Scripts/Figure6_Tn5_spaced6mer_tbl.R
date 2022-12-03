library(parallel)
library(bigWig)
source("https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_hcsearch.R")
#Load in all spaced 3-mer masks
masks = read.table('all_possible_spaced_3mers_46bp.txt', skip = 1)
#Assign each mask to a list object in nextlst
nextlst = vector(mode = 'list', length = nrow(masks))
for (i in 1:nrow(masks)) {
  nextlst[[i]][1] = masks[i,1]
}
#Arguments for seqOutBias
seqOutBias.args = "hg38.fa C1_gDNA_rep1.bam --read-size=76 --strand-specific --custom-shift=4,-4"
#Command to call seqOutBias
sqcmd = "seqOutBias"
#Prefix for output
prefix = "Tn5_"
#Run seqOutBias for each mask and get their .tbl files
scores = mclapply(nextlst, function(cutmask) {
  bw.paths = run.cutmask(cutmask, seqOutBias.args, sqcmd=sqcmd, prefix=prefix, cleanup = FALSE)
}, mc.cores = 1)