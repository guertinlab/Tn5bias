library(parallel)
library(bigWig)
source("https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_hcsearch.R")
start_mask = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
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
      vi[k+5] = 1 # 6th N
      result[[n]] <- vec.to.mask(vi)
    }
  }
  if (n == 0) return(NULL)
  #Remove last line so you don't get an extra N space (mappable bases = mer-1)
  zz <- (n-5)
  result[1:zz]
}
nextlst = neighbors(mask.to.vec(start_mask))
seqOutBias.args = "hg38.fa C1_gDNA_rep1.bam --read-size=76 --strand-specific --custom-shift=4,-4"
sqcmd = "seqOutBias"
prefix = "Tn5_"
scores = mclapply(nextlst, function(cutmask) {
  bw.paths = run.cutmask(cutmask, seqOutBias.args, sqcmd=sqcmd, prefix=prefix, cleanup = FALSE)
}, mc.cores = 1)