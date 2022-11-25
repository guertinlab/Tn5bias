library(parallel)
library(bigWig)
source("https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_hcsearch.R")

system('mkdir EnzymeBias_masks')
setwd('EnzymeBias_masks')
neighbors <- function(vecmask) {
  # generate list of neighboring masks that differ
  # by the addition of a single unmasked position
  result <- vector(mode="list", length=length(vecmask))
  n <- 0
  for (k in 1:length(vecmask)) {
    if (vecmask[k] == -1) { # X
      n <- n + 1
      vi = vecmask
      vi[k] = 1 # N
      vi[k+1] = 1 # 2nd N
      vi[k+2] = 1 # 3rd N
      vi[k+3] = 1 # 4th N
      vi[k+4] = 1 # 5th N
      result[[n]] <- vec.to.mask(vi)
    }
  }
  if (n == 0) return(NULL)
  #Remove last line so you don't get an extra N space (mappable bases = mer-1)
  zz <- (n-4)
  result[1:zz]
}

###Make Tn5 masks-we need a wider range because we will need to shift these by 4 (this will take hours)
start_mask = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
nextlst = neighbors(mask.to.vec(start_mask))
seqOutBias.args = "../hg38.fa ../C1_gDNA_rep1.bam --read-size=76 --strand-specific --custom-shift=4,-4"
sqcmd = "seqOutBias"
prefix = "Tn5_"

masks = mclapply(nextlst, function(cutmask) {
  bw.paths = run.cutmask(cutmask, seqOutBias.args, sqcmd=sqcmd, prefix=prefix, cleanup = FALSE)
}, mc.cores = 3)

###Make DNase masks (this will take hours)
start_mask = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
nextlst = neighbors(mask.to.vec(start_mask))
seqOutBias.args = "../hg38.fa ../DNase_Naked.bam --read-size=76 --strand-specific"
sqcmd = "seqOutBias"
prefix = "DNase_"

masks = mclapply(nextlst, function(cutmask) {
  bw.paths = run.cutmask(cutmask, seqOutBias.args, sqcmd=sqcmd, prefix=prefix, cleanup = FALSE)
}, mc.cores = 3)

###Make Benzonase masks (this will take hours)
start_mask = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
nextlst = neighbors(mask.to.vec(start_mask))
seqOutBias.args = "../mm39.fa ../mm39_liver_Benzonase.bam --read-size=76 --strand-specific"
sqcmd = "seqOutBias"
prefix = "Benzonase_"

masks = mclapply(nextlst, function(cutmask) {
  bw.paths = run.cutmask(cutmask, seqOutBias.args, sqcmd=sqcmd, prefix=prefix, cleanup = FALSE)
}, mc.cores = 3)

###Make Cyanase masks (this will take hours)
start_mask = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
nextlst = neighbors(mask.to.vec(start_mask))
seqOutBias.args = "../mm39.fa ../mm39_liver_Cyanase.bam --read-size=76 --strand-specific"
sqcmd = "seqOutBias"
prefix = "Cyanase_"

masks = mclapply(nextlst, function(cutmask) {
  bw.paths = run.cutmask(cutmask, seqOutBias.args, sqcmd=sqcmd, prefix=prefix, cleanup = FALSE)
}, mc.cores = 3)

###Make MNase masks (this will take hours)
start_mask = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
nextlst = neighbors(mask.to.vec(start_mask))
seqOutBias.args = "../mm39.fa ../mm39_liver_MNase.bam --read-size=76 --strand-specific"
sqcmd = "seqOutBias"
prefix = "MNase_"

masks = mclapply(nextlst, function(cutmask) {
  bw.paths = run.cutmask(cutmask, seqOutBias.args, sqcmd=sqcmd, prefix=prefix, cleanup = FALSE)
}, mc.cores = 3)