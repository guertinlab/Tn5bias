library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
#Load in hg38 plus sequences from Figure 5
load('../Figure5/hg38_Plus_motif_sequences.Rdata')
#Break sequences into 7bp chunks (for 7-mer)
TFseq_plus = lapply(TFseq_plus, mer.positions, mermask = "NNNNNNN")
#Find frequency of each possible 7-mer for each position
TFseq_plus = lapply(TFseq_plus, position.frequencies, mermask = "NNNNNNN")
save(TFseq_plus, file = 'Tn5_TF_position_frequencies_7mer_plus.Rdata')

#Determine minus strand 7-mer frequencies
load('../Figure5/hg38_Minus_motif_sequences.Rdata')
#Break sequences into 7bp chunks (for 7-mer)
TFseq_minus = lapply(TFseq_minus, mer.positions, mermask = "NNNNNNN")
#Find frequency of each possible 7-mer for each position
TFseq_minus = lapply(TFseq_minus, position.frequencies, mermask = "NNNNNNN")
save(TFseq_minus, file = 'Tn5_TF_position_frequencies_7mer_minus.Rdata')
