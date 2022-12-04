library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
system('mkdir spaced_3mer')
setwd('spaced_3mer')
system('mkdir plus')
setwd('plus')
#Determine the unique spacings of 3mers in 46bp
factorfiles <- read.table('../../Tn5_kmercount/all_possible_spaced_6mers_46bp.txt', skip = 1)
factorfiles = factorfiles[,1]
unique_factorfiles_spacing = factorfiles
unique_factorfiles_spacing = strsplit(unique_factorfiles_spacing, 'NNN')
unique_factorfiles = NULL
for (i in 1:length(unique_factorfiles_spacing)) {
  unique_factorfiles[i] = unique_factorfiles_spacing[[i]][2]
}
unique_factorfiles_spacing = gsub('^X', 'NNNX', unique_factorfiles)
unique_factorfiles_spacing = gsub('X$', 'XNNN', unique_factorfiles_spacing)
unique_factorfiles_spacing = unique(unique_factorfiles_spacing)
rm(unique_factorfiles)
###############################################################
load('../../../Figure5/hg38_Plus_motif_sequences.Rdata')
for (i in 1:length(unique_factorfiles_spacing)) {
  TFseq_plus_s = lapply(TFseq_plus, mer.positions,
                        mermask = unique_factorfiles_spacing[i])
  TFseq_plus_s = lapply(TFseq_plus_s, position.frequencies,
                        mermask = unique_factorfiles_spacing[i])
  save(TFseq_plus_s, file = paste(unique_factorfiles_spacing[i],
                                  '_TFseq_plus_posfreqs.Rdata', sep = ''))
}
###############################################################
system('mkdir ../minus')
setwd('../minus')
load('../../../Figure5/hg38_Minus_motif_sequences.Rdata')
for (i in 1:length(unique_factorfiles_spacing)) {
  TFseq_minus_s = lapply(TFseq_minus, mer.positions,
                         mermask = unique_factorfiles_spacing[i])
  TFseq_minus_s = lapply(TFseq_minus_s, position.frequencies,
                         mermask = unique_factorfiles_spacing[i])
  save(TFseq_minus_s, file = paste(unique_factorfiles_spacing[i],
                                   '_TFseq_minus_posfreqs.Rdata', sep = ''))
}
setwd('../../')
