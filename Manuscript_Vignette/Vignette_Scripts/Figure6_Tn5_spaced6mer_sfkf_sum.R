library(data.table)
library('parallel')
library('foreach')
library('doParallel')
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
registerDoParallel(20)
setwd('../sfkf_sum')
options(scipen = 100)

RDSlist = list.files('../plus_HPC_output')
TF_plus_scalefactors_kmerfreq <- vector(mode = "list", length = length(RDSlist))
TF_plus_scalefactors_kmerfreq = foreach (i = 1:length(RDSlist)) %dopar% {
  RDSstore = readRDS(paste( '../plus_HPC_output/', RDSlist[i], sep = ''))
  RDSstore = scaledkmerfreq.sum.plus(RDSstore)
  return(RDSstore)}

names(TF_plus_scalefactors_kmerfreq) = RDSlist

save(TF_plus_scalefactors_kmerfreq, 
     file = 'Tn5_TF_plus_pre_input_spaced6mers.Rdata')


RDSlist = list.files('../minus_HPC_output')
TF_minus_scalefactors_kmerfreq <- vector(mode = "list", length = length(RDSlist))
TF_minus_scalefactors_kmerfreq = foreach (i = 1:length(RDSlist)) %dopar% {
  RDSstore = readRDS(paste( '../minus_HPC_output/', RDSlist[i], sep = ''))
  RDSstore = scaledkmerfreq.sum.minus(RDSstore)
  return(RDSstore)}

names(TF_minus_scalefactors_kmerfreq) = RDSlist

save(TF_minus_scalefactors_kmerfreq, 
     file = 'Tn5_TF_minus_pre_input_spaced6mers.Rdata')
