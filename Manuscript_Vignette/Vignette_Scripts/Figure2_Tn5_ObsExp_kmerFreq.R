source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
library(data.table)
library(lattice)
library(stringr)
options(scipen=10000)
#
#
uppercasenames <- list('../Figure1/C1_gDNA_rep1_plus.fasta', '../Figure1/C1_gDNA_rep1_minus_RC.fasta')
uplist <- lapply(uppercasenames, uppercase)  
plus_Tn5 = as.data.frame(uplist[1])

minus_Tn5 = as.data.frame(uplist[2])
rm(uplist, uppercasenames)
#
#######################################################################################################
#Make df of positions 10-12 (-6:-4 from cutsite)
plus_Tn5_pos10_12 = data.table(str_split_fixed(plus_Tn5[,1], '', 13))
plus_Tn5_pos10_12 = plus_Tn5_pos10_12[,10:12]
plus_Tn5_pos10_12[, kmer := paste0(V10, V11, V12)]
plus_Tn5_pos10_12 = plus_Tn5_pos10_12$kmer
rm(plus_Tn5)

#Make df of positions 10-12 (-6:-4 from cutsite)
minus_Tn5_pos10_12 = data.table(str_split_fixed(minus_Tn5[,1], '', 13))
minus_Tn5_pos10_12 = minus_Tn5_pos10_12[,10:12]
minus_Tn5_pos10_12[, kmer := paste0(V10, V11, V12)]
minus_Tn5_pos10_12 = minus_Tn5_pos10_12$kmer
rm(minus_Tn5)

#Make all possible 3mers:
mertable = expand.grid(rep(list(c('A','C','G','T')), 3))
mertable = data.frame(apply(mertable, 1 , paste, collapse = ""))

#Count number of each 3mer in each position:
for (i in 1:length(mertable[,1])) {
  mertable[i,2] = length(grep(mertable[i,1], plus_Tn5_pos10_12))
  mertable[i,3] = length(grep(mertable[i,1], minus_Tn5_pos10_12))
}

colnames(mertable) = c('kmer', 'upPlus', 'upMinus')

save(mertable, file = 'mertable_kmer_counts.Rdata')
load('mertable_kmer_counts.Rdata')

mertable[,4] = mertable[,2] + mertable[,3]
colnames(mertable)[4] = c('upTotal')
mertable[,4] = mertable[,4]/sum(mertable[,4])
##Load in tn5 bias transfac
transfac.tn5 = read.table('../Figure1/Tn5_sepcat_bias.transfac', skip = 3, fill = TRUE)
transfac.tn5 = transfac.tn5[-c(32), -c(1)]
colnames(transfac.tn5) = c('A', 'C', 'G', 'T')
transfac.tn5 = t(apply(transfac.tn5, 1, function(x) x / sum(x)))
###Calculate expected upstream values based on nucleotide prevalence:
expected_upstream_nuc_prevalence = expand.grid(rep(list(c('A','C','G','T')), 3))

for (i in 1:nrow(expected_upstream_nuc_prevalence)) {
  expected_upstream_nuc_prevalence[i,4] = 
    transfac.tn5[10 ,which(colnames(transfac.tn5) == expected_upstream_nuc_prevalence[i,1])]
  expected_upstream_nuc_prevalence[i,5] = 
    transfac.tn5[11 ,which(colnames(transfac.tn5) == expected_upstream_nuc_prevalence[i,2])]
  expected_upstream_nuc_prevalence[i,6] = 
    transfac.tn5[12 ,which(colnames(transfac.tn5) == expected_upstream_nuc_prevalence[i,3])]
}

#Multiply these together to get each 3-mers expected prevalence
expected_upstream_nuc_prevalence[,7] = expected_upstream_nuc_prevalence[,4] *
  expected_upstream_nuc_prevalence[,5] * 
  expected_upstream_nuc_prevalence[,6]

#calculate upstream observed / expected:
mertable$upOE = mertable$upTotal / expected_upstream_nuc_prevalence$V7
#Correct the observed divided by expected k-mer counts by the 50bp upstream counts
load('Tn5_100bp_observedexpected_kmercounts.Rdata')
mertable$CorrupOE = mertable$upOE*(1/Tn5_kmer_count_correction$OE50bp)
mertable$CorrupOE = log2(mertable$CorrupOE)
mertable$kmer = factor(mertable$kmer, levels = mertable$kmer)
save(mertable, file = 'Figure2A_Tn5_upstream_OE_kmers.Rdata')
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
