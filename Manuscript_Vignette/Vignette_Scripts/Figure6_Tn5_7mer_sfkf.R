library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
library('parallel')
library('doParallel')
options(scipen = 100)

#Import kmer counts and compute scale factors for each....
factorfiles <- list.files('Tn5_kmercount')
factorfiles = factorfiles[grep('_scale_factors.txt', factorfiles)]
gsub_factorfiles = gsub('C', '', factorfiles)
factorfiles = factorfiles[grep('NNNNNNN', gsub_factorfiles)]
gsub_factorfiles = gsub_factorfiles[grep('NNNNNNN', gsub_factorfiles)]
rm(gsub_factorfiles)
scalefactors <- vector(mode = "list", length = length(factorfiles))
names(scalefactors) <- factorfiles

for (i in 1:length(factorfiles)) {
  scalefactors[[i]] <- scalefactor.func(paste('Tn5_kmercount/',
                                              factorfiles[i], sep = ''))
}

##7mer Plus scale factor by k-mer frequency multiplication
load('Tn5_TF_position_frequencies_7mer_plus.Rdata')

registerDoParallel(3)
TF_scalefactors_kmerfreq <- vector(mode = "list", length = length(TFseq_plus))
names(TF_scalefactors_kmerfreq) = names(TFseq_plus)

foreach (i = 1:length(TF_scalefactors_kmerfreq)) %dopar% 
  {TF_scalefactors_kmerfreq[[i]] = scalefactor.by.kmerfrequency(scalefactors = scalefactors,
                                                                kmerfrequency = TFseq_plus[[i]], 
                                                                tfname = names(TFseq_plus[i]))}

TF_sfkf_7mer_minus <- vector(mode = "list", length = length(TFseq_minus))
names(TF_sfkf_7mer_minus) = names(TFseq_minus)
rm(TFseq_minus)

RDSlist = list.files('./')
RDSlist = RDSlist[grep('.rds', RDSlist)]
RDSlist = RDSlist[grep('plus', RDSlist)]
RDSlist = RDSlist[grep(paste(substr(factorfiles,5,nchar(factorfiles)-18),
                             collapse="|"),RDSlist)]

for (i in 1:length(RDSlist)) {
  TF_sfkf_7mer_plus[[which(names(TF_sfkf_7mer_plus) ==
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 56))]][[grep(substr(RDSlist[i],
                    nchar(RDSlist[i]) - 50, nchar(RDSlist[i])- 4), factorfiles)]] =
                    readRDS(paste(RDSlist[i], sep = ''))
  
  names(TF_sfkf_7mer_plus[[which(names(TF_sfkf_7mer_plus) ==
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 56))]])[grep(substr(RDSlist[i],
                    nchar(RDSlist[i]) - 50, nchar(RDSlist[i])- 4), factorfiles)] =
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 4)
  
  TF_sfkf_7mer_plus[[which(names(TF_sfkf_7mer_plus) ==
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 56))]][[grep(substr(RDSlist[i],
                    nchar(RDSlist[i]) - 50, nchar(RDSlist[i])- 4), factorfiles)]] =
                    scaledkmerfreq.sum.plus(TF_sfkf_7mer_plus[[which(names(TF_sfkf_7mer_plus) ==
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 56))]][[grep(substr(RDSlist[i],
                    nchar(RDSlist[i]) - 50, nchar(RDSlist[i])- 4), factorfiles)]])
}

#Shift the values so that all predictions are aligned
input_correction = data.frame(matrix(nrow = 40, ncol = 2))
input_correction[1:40,1] = 1:40
input_correction[1:40,2] = 156:195

plus_pred_input = vector(mode = 'list', length = length(TF_sfkf_7mer_plus))
names(plus_pred_input) <- names(TF_sfkf_7mer_plus)
for (i in 1:length(TF_sfkf_7mer_plus)) {
  for (p in 1:length(TF_sfkf_7mer_plus[[i]])) {
    plus_pred_input[[i]][[p]] = TF_sfkf_7mer_plus[[i]][[p]][input_correction[p,1]:input_correction[p,2]]
  }
  plus_pred_input[[i]] = data.frame(plus_pred_input[[i]])
  names(plus_pred_input[[i]]) = substr(names(TF_sfkf_7mer_plus[[i]]),
  nchar(names(TF_sfkf_7mer_plus[[i]])) - 46, nchar(names(TF_sfkf_7mer_plus[[i]])))
}

###Add x axis to input for future alignment
for (i in 1:length(plus_pred_input)) {
  for (p in 1:length(plus_pred_input[[i]])) {
    plus_pred_input[[i]][[p]] = as.data.frame(plus_pred_input[[i]][[p]])
    plus_pred_input[[i]][[p]][,2] = seq(-73,80,1)
    colnames(plus_pred_input[[i]][[p]]) = c('sfkf', 'x')
  }
}

save(plus_pred_input, file = 'Tn5_TF_plus_pre_input_7mers.Rdata')

################################################################################
##7mer minus scale factor by k-mer frequency multiplication
load('Tn5_TF_position_frequencies_7mer_minus.Rdata')


registerDoParallel(3)
TF_scalefactors_kmerfreq <- vector(mode = "list", length = length(TFseq_minus))
names(TF_scalefactors_kmerfreq) = names(TFseq_minus)

foreach (i = 1:length(TF_scalefactors_kmerfreq)) %dopar% 
  {TF_scalefactors_kmerfreq[[i]] = scalefactor.by.kmerfrequency(scalefactors = scalefactors,
                                                                kmerfrequency = TFseq_minus[[i]],
                                                                tfname = names(TFseq_minus[i]))}


TF_sfkf_7mer_minus <- vector(mode = "list", length = length(TFseq_minus))
names(TF_sfkf_7mer_minus) = names(TFseq_minus)
rm(TFseq_minus)

RDSlist = list.files('./')
RDSlist = RDSlist[grep('.rds', RDSlist)]
RDSlist = RDSlist[grep('minus', RDSlist)]
RDSlist = RDSlist[grep(paste(substr(factorfiles,5,nchar(factorfiles)-18),
                             collapse="|"),RDSlist)]

for (i in 1:length(RDSlist)) {
  TF_sfkf_7mer_minus[[which(names(TF_sfkf_7mer_minus) ==
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 56))]][[grep(substr(RDSlist[i],
                    nchar(RDSlist[i]) - 50, nchar(RDSlist[i])- 4), factorfiles)]] =
                    readRDS(paste(RDSlist[i], sep = ''))
  
  names(TF_sfkf_7mer_minus[[which(names(TF_sfkf_7mer_minus) ==
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 56))]])[grep(substr(RDSlist[i],
                    nchar(RDSlist[i]) - 50, nchar(RDSlist[i])- 4), factorfiles)] =
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 4)
  
  TF_sfkf_7mer_minus[[which(names(TF_sfkf_7mer_minus) ==
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 56))]][[grep(substr(RDSlist[i],
                    nchar(RDSlist[i]) - 50, nchar(RDSlist[i])- 4), factorfiles)]] =
                    scaledkmerfreq.sum.minus(TF_sfkf_7mer_minus[[which(names(TF_sfkf_7mer_minus) ==
                    substr(RDSlist[i], 1, nchar(RDSlist[i]) - 56))]][[grep(substr(RDSlist[i],
                    nchar(RDSlist[i]) - 50, nchar(RDSlist[i])- 4), factorfiles)]])
}

#Shift the values so that all predictions are aligned
input_correction = data.frame(matrix(nrow = 40, ncol = 2))
input_correction[1:40,1] = 1:40
input_correction[1:40,2] = 156:195

minus_pred_input = vector(mode = 'list', length = length(TF_sfkf_7mer_minus))
names(minus_pred_input) <- names(TF_sfkf_7mer_minus)
for (i in 1:length(TF_sfkf_7mer_minus)) {
  for (p in 1:length(TF_sfkf_7mer_minus[[i]])) {
    minus_pred_input[[i]][[p]] = TF_sfkf_7mer_minus[[i]][[p]][input_correction[p,1]:input_correction[p,2]]
  }
  minus_pred_input[[i]] = data.frame(minus_pred_input[[i]])
  names(minus_pred_input[[i]]) = substr(names(TF_sfkf_7mer_minus[[i]]),
        nchar(names(TF_sfkf_7mer_minus[[i]])) - 46, nchar(names(TF_sfkf_7mer_minus[[i]])))
}

###Add x axis to input for future alignment
for (i in 1:length(minus_pred_input)) {
  for (p in 1:length(minus_pred_input[[i]])) {
    minus_pred_input[[i]][[p]] = as.data.frame(minus_pred_input[[i]][[p]])
    minus_pred_input[[i]][[p]][,2] = seq(-73,80,1)
    colnames(minus_pred_input[[i]][[p]]) = c('sfkf', 'x')
  }
}

save(minus_pred_input, file = 'Tn5_TF_minus_pre_input_7mers.Rdata')
