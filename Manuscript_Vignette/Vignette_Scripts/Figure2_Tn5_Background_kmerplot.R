source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
library(data.table)
library(ggplot2)
options(scipen = 100)
#Convert all sequences to uppercase
uppercasenames <- list('C1_gDNA_rep1_100bp_sepcat.fasta', 'C1_gDNA_rep1_50bp_sepcat.fasta')
uplist <- lapply(uppercasenames, uppercase)  
Tn5_100bp_sepcat = as.data.frame(uplist[1])
Tn5_50bp_sepcat = as.data.frame(uplist[2])
#Determine pswm of each:
Tn5_100bp_sepcat.pswm = transfac.func.2(Tn5_100bp_sepcat[,1], 3)
Tn5_100bp_sepcat.pswm = t(apply(Tn5_100bp_sepcat.pswm, 1, function(x) x / sum(x)))

Tn5_50bp_sepcat.pswm = transfac.func.2(Tn5_50bp_sepcat[,1], 3)
Tn5_50bp_sepcat.pswm = t(apply(Tn5_50bp_sepcat.pswm, 1, function(x) x / sum(x)))
#Make df of positions 50-47bp upstream (-50:-47 from cutsite)
Tn5_pos100_98 = data.table(substr(uplist[[1]][,1], 1, 3))
Tn5_pos100_98 = Tn5_pos100_98$V1
Tn5_pos50_48 = data.table(substr(uplist[[2]][,1], 1, 3))
Tn5_pos50_48 = Tn5_pos50_48$V1
#Make all possible 3mers:
mertable = expand.grid(rep(list(c('A','C','G','T')), 3))
mertable = data.frame(apply(mertable, 1 , paste, collapse = ""))
#Count number of each 3mer in each position:
for (i in 1:length(mertable[,1])) {
  mertable[i,2] = length(grep(mertable[i,1], Tn5_pos100_98))
  mertable[i,3] = length(grep(mertable[i,1], Tn5_pos50_48))
}
#Calculate percent frequency of each k-mer for each position
mertable[,4] = mertable[,2]/sum(mertable[,2])
mertable[,5] = mertable[,3]/sum(mertable[,3])
#Change column names
colnames(mertable) = c('kmer', 'count100bp', 'count50bp', 'percent100bp', 'percent50bp')
#Calculate expected k-mer frequency for 100bp, based on the PSWM created above
pswm_100bp_Tn5 = Tn5_100bp_sepcat.pswm
colnames(pswm_100bp_Tn5) = c('A', 'C', 'G', 'T')
###Calculate expected upstream values based on nucleotide prevalence:
expected_100bp_nuc_prevalence = expand.grid(rep(list(c('A','C','G','T')), 3))
for (i in 1:nrow(expected_100bp_nuc_prevalence)) {
  expected_100bp_nuc_prevalence[i,4] = 
    pswm_100bp_Tn5[1 ,which(colnames(pswm_100bp_Tn5) == expected_100bp_nuc_prevalence[i,1])]
  expected_100bp_nuc_prevalence[i,5] = 
    pswm_100bp_Tn5[2 ,which(colnames(pswm_100bp_Tn5) == expected_100bp_nuc_prevalence[i,2])]
  expected_100bp_nuc_prevalence[i,6] = 
    pswm_100bp_Tn5[3 ,which(colnames(pswm_100bp_Tn5) == expected_100bp_nuc_prevalence[i,3])]
}
expected_100bp_nuc_prevalence[,7] = 
  expected_100bp_nuc_prevalence[,4] * expected_100bp_nuc_prevalence[,5] * expected_100bp_nuc_prevalence[,6]
#calculate 100bp observed / expected:
mertable$OE100bp = mertable$percent100bp / expected_100bp_nuc_prevalence$V7
#Repeat for 50bp upstream
#Calculate expected k-mer frequency for 50bp, based on the PSWM created above
pswm_50bp_Tn5 = Tn5_50bp_sepcat.pswm
colnames(pswm_50bp_Tn5) = c('A', 'C', 'G', 'T')
###Calculate expected upstream values based on nucleotide prevalence:
expected_50bp_nuc_prevalence = expand.grid(rep(list(c('A','C','G','T')), 3))
for (i in 1:nrow(expected_50bp_nuc_prevalence)) {
  expected_50bp_nuc_prevalence[i,4] = 
    pswm_50bp_Tn5[1 ,which(colnames(pswm_50bp_Tn5) == expected_50bp_nuc_prevalence[i,1])]
  expected_50bp_nuc_prevalence[i,5] = 
    pswm_50bp_Tn5[2 ,which(colnames(pswm_50bp_Tn5) == expected_50bp_nuc_prevalence[i,2])]
  expected_50bp_nuc_prevalence[i,6] = 
    pswm_50bp_Tn5[3 ,which(colnames(pswm_50bp_Tn5) == expected_50bp_nuc_prevalence[i,3])]
}
expected_50bp_nuc_prevalence[,7] = 
  expected_50bp_nuc_prevalence[,4] * expected_50bp_nuc_prevalence[,5] * expected_50bp_nuc_prevalence[,6]
#calculate 50bp observed / expected:
mertable$OE50bp = mertable$percent50bp / expected_50bp_nuc_prevalence$V7
#Calculate 1/50bp * 100bp Observed/Expected values 
mertable$Inv50OE100bp = mertable$OE100bp*(1/mertable$OE50bp)
#Plot this comparison
pdf(file = 'Tn5_100bpCorrectedby50bp_kmers.pdf', width=12, height=9)
ggplot(data=mertable, aes(x=kmer, y=Inv50OE100bp)) +
  geom_bar(stat="identity", fill = 'Black') +
  ylab("Ratio observed / expected") + xlab("") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 16, colour = 'black'),
        axis.text.y = element_text(size = 16, color = 'black'),
        axis.title.y = element_text(size = 24)) +
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=1)
dev.off()
#Save the mertable df to use as a correction for observed/expected k-mer frequencies
Tn5_kmer_count_correction = mertable
save(Tn5_kmer_count_correction, file='Tn5_100bp_observedexpected_kmercounts.Rdata')
