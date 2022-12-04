library(data.table)
library(doParallel)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
train = c( 'ASCL1', 'ATF3', 'CLOCK', 'CTCF',
           'DLX2', 'EBF1', 'EGR1', 'ELF1',
           'FOS', 'FOXA1', 'FOXK2', 'GATA2',
           'GLIS1', 'IRF1', 'JUN', 'LEF1',
           'MEIS2', 'MLX', 'MYC', 'PPARG', 
           'RUNX1', 'SPIB', 'SRY', 'STAT1', 'TGIF1')
#Import kmer counts and compute scale factors for each
factorfiles = list.files('../Figure1')
factorfiles = factorfiles[c(grep('DNase', factorfiles))]
factorfiles = factorfiles[c(grep('_scale_factors.txt', factorfiles))]
factorfiles = factorfiles[grep('NNNNNXXXXXCXXXXXXXXXX', factorfiles):
                            grep('XXXXXXXXXXCXXXXXNNNNN', factorfiles)]
scalefactors = vector(mode = "list", length = length(factorfiles))
names(scalefactors) = factorfiles
for (i in 1:length(factorfiles)) {
  scalefactors[[i]] <- scalefactor.func(paste('../Figure1/', factorfiles[i], sep = ''))
}


#Import TF FIMO coordinates to extract the sequences from the FASTA
FIMOfiles = list.files('./')
FIMOfiles = FIMOfiles[grep('plus_400k_fimo', FIMOfiles)]

FIMOfiles = FIMOfiles[which(substr(FIMOfiles,1,nchar(FIMOfiles)-22) %in% train)]
TFmotifs = vector(mode = "list", length = length(FIMOfiles))
names(TFmotifs) = substr(FIMOfiles, 1, (nchar(FIMOfiles)-14))

#Import FIMO files for genomic coordinates and expand a 'window' around these sites (using a central base)
for (i in 1:length(TFmotifs)) {
  TFmotifs[[i]] <- FIMO.to.BED(paste('./', FIMOfiles[i], sep = ''), window = 100)
}

#Ensure that no start positions are off the chromosome
for (i in 1:length(TFmotifs)) {
  print(TFmotifs[[i]][which.min(TFmotifs[[i]]$start),])
}

#Write these coordinates into a bed file for use with 'getfasta' from bedtools
for (i in 1:length(TFmotifs)) {
  write.table(TFmotifs[[i]], file = paste(names(TFmotifs[i]), '_hg38coords.bed', sep = ''),
              row.names=FALSE, sep="\t", quote = FALSE, col.names = FALSE)
}

#Run getfasta on genomic coordinates to get sequences for window
for (i in 1:length(TFmotifs)) {
  system(paste('bedtools getfasta -bedOut -s -fi ../Figure1/hg38.fa -bed ', 
               names(TFmotifs[i]), '_hg38coords.bed > ',
               names(TFmotifs[i]), '_hg38seq.bed', sep = ''))
}
#load in bed files with sequences
TFseq_plus <- vector(mode = "list", length = length(TFmotifs))
names(TFseq_plus) <- names(TFmotifs)
for (i in 1:length(TFmotifs)) {
  TFseq_plus[[i]] <- fread(paste(names(TFmotifs[i]), '_hg38seq.bed', sep = ''))
}
#Make sure you still have 400k rows
for (i in 1:length(TFseq_plus)) {
  print(nrow(TFseq_plus[[i]]))
}
#Make all sequences upper case
TFseq_plus = lapply(TFseq_plus, function(x) {toupper(x$V7)})
save(TFseq_plus, file = 'hg38_Plus_motif_sequences.Rdata')
#Break sequences into 5bp chunks (for 5mer)
TFseq_plus = lapply(TFseq_plus, mer.positions, mermask = "NNNNN")
#Find frequency of each possible 5mer for each position
TFseq_plus = lapply(TFseq_plus, position.frequencies, mermask = 'NNNNN')
#Multiply 1/scale_factor by its corresponding kmer frequency for each TF motif. This will produce an
# .rds file for each TF/mask combo (scale factor * k-mer frequency - sfkf)
registerDoParallel(3)
TF_scalefactors_kmerfreq <- vector(mode = "list", length = length(TFseq_plus))
names(TF_scalefactors_kmerfreq) = names(TFseq_plus)

foreach (i = 1:length(TF_scalefactors_kmerfreq)) %dopar% {TF_scalefactors_kmerfreq[[i]] = 
                               scalefactor.by.kmerfrequency(scalefactors = scalefactors,
                               kmerfrequency = TFseq_plus[[i]], tfname = names(TFseq_plus[i]))}

#Next we load in each sfkf combo to the same list
TF_sfkf_5mer_plus <- vector(mode = "list", length = length(TFseq_plus))
names(TF_sfkf_5mer_plus) = names(TFseq_plus)

RDSlist = list.files('./')
RDSlist = RDSlist[grep('plus', RDSlist)]
RDSlist = RDSlist[grep('.rds', RDSlist)]

for (i in 1:length(RDSlist)) {
  TF_sfkf_5mer_plus[[which(names(TF_sfkf_5mer_plus) == 
          substr(RDSlist[i], 1, nchar(RDSlist[i]) - 48))]][[grep(substr(RDSlist[i],
          nchar(RDSlist[i]) - 40, nchar(RDSlist[i])- 4),
          factorfiles)]] = readRDS(paste( './', RDSlist[i], sep = ''))
          names(TF_sfkf_5mer_plus[[which(names(TF_sfkf_5mer_plus) == 
          substr(RDSlist[i], 1, nchar(RDSlist[i]) - 48))]])[grep(substr(RDSlist[i],
          nchar(RDSlist[i]) - 40, nchar(RDSlist[i])- 4),
          factorfiles)] =  substr(RDSlist[i], 1, nchar(RDSlist[i]) - 4)
          TF_sfkf_5mer_plus[[which(names(TF_sfkf_5mer_plus) == 
          substr(RDSlist[i], 1, nchar(RDSlist[i]) - 48))]][[grep(substr(RDSlist[i],
          nchar(RDSlist[i]) - 40, nchar(RDSlist[i])- 4),
          factorfiles)]] = scaledkmerfreq.sum.plus(TF_sfkf_5mer_plus[[which(names(TF_sfkf_5mer_plus)
          == substr(RDSlist[i], 1,
          nchar(RDSlist[i]) - 48))]][[grep(substr(RDSlist[i], nchar(RDSlist[i]) - 40,
                                                                                                                                                                                                        nchar(RDSlist[i])- 4), factorfiles)]])
}
Plus_pre_input = TF_sfkf_5mer_plus

###############################################################################
#Load in unscaled composites
load('DNase_plus_unscaled_compositelist.Rdata')
DNase_unscaled_composites = do.call(rbind, DNase_plus_unscaled_compositelist)
#Trim unscaled composites to the size of predicted area
DNase_unscaled_composites = DNase_unscaled_composites[-c(which(DNase_unscaled_composites$x < -90.5 |
                                                                 DNase_unscaled_composites$x > 90.5)),]
#Normalize unscaled composites to the same value for all composites
DNase_unscaled_composites$est = DNase_unscaled_composites$est/mean(DNase_unscaled_composites$correction)
#Subset out training motifs
DNase_unscaled_composites = DNase_unscaled_composites[which(substr(DNase_unscaled_composites$factor,1,
                            nchar(DNase_unscaled_composites$factor)-8) %in% train),]
#Shift pre predictor variables. Each masks' position shifts by the mask's offset
input_correction = data.frame(matrix(nrow = 16, ncol = 2))
input_correction[1:16,1] = 1:16
input_correction[1:16,2] = 182:197

Plus_pred_input = vector(mode = 'list', length = length(Plus_pre_input))
names(Plus_pred_input) <- names(Plus_pre_input)
for (i in 1:length(Plus_pre_input)) {
  for (p in 1:length(Plus_pre_input[[i]])) {
    Plus_pred_input[[i]][[p]] = Plus_pre_input[[i]][[p]][input_correction[p,1]:input_correction[p,2]]
  }
  Plus_pred_input[[i]] = data.frame(Plus_pred_input[[i]])
}

#Assign names
for (i in 1:length(Plus_pred_input)) {
  names(Plus_pred_input[[i]]) = names(Plus_pre_input[[i]])  
}
#Make all mask names the same:
for (i in 1:length(Plus_pred_input)) {
  names(Plus_pred_input[[i]]) = substr(names(Plus_pred_input[[i]]),
                    nchar(names(Plus_pred_input[[i]]))-36, nchar(names(Plus_pred_input[[i]])))
}
#Turn plus_pred_input into DF
Plus_pred_input = cbind(seq(-90.5,90.5,1), DNase_unscaled_composites$est, do.call(rbind, Plus_pred_input))
colnames(Plus_pred_input)[1:2] = c('X', 'unscaled')
Plus_pred_input = Plus_pred_input[which(Plus_pred_input$X >= -39.5 & Plus_pred_input$X <= 39.5),]
#This is the plus input for the rule ensemble
save(Plus_pred_input, file = 'DNase_TF_plus_pre_input.Rdata')

###############################################################################
#Repeat for Minus strand motifs
#Import TF FIMO coordinates to extract the sequences from the FASTA
FIMOfiles = list.files('./')
FIMOfiles = FIMOfiles[grep('rm_minus_400k_fimo.txt', FIMOfiles)]
FIMOfiles = FIMOfiles[which(substr(FIMOfiles,1,nchar(FIMOfiles)-23) %in% train)]
TFmotifs = vector(mode = "list", length = length(FIMOfiles))
names(TFmotifs) = substr(FIMOfiles, 1, (nchar(FIMOfiles)-14))

#Import FIMO files for genomic coordinates and expand a 'window' around these sites (using a central base)
for (i in 1:length(TFmotifs)) {
  TFmotifs[[i]] <- FIMO.to.BED(paste('./', FIMOfiles[i], sep = ''), window = 100)
}

#Ensure that no start positions are off the chromosome
for (i in 1:length(TFmotifs)) {
  print(TFmotifs[[i]][which.min(TFmotifs[[i]]$start),])
}

#Write these coordinates into a bed file for use with 'getfasta' from bedtools
for (i in 1:length(TFmotifs)) {
  write.table(TFmotifs[[i]], file = paste(names(TFmotifs[i]), '_hg38coords.bed', sep = ''),
              row.names=FALSE, sep="\t", quote = FALSE, col.names = FALSE)
}

#Run getfasta on genomic coordinates to get sequences for window
for (i in 1:length(TFmotifs)) {
  system(paste('bedtools getfasta -bedOut -s -fi ../Figure1/hg38.fa -bed ', 
               names(TFmotifs[i]), '_hg38coords.bed > ',
               names(TFmotifs[i]), '_hg38seq.bed', sep = ''))
}
#load in bed files with sequences
TFseq_minus <- vector(mode = "list", length = length(TFmotifs))
names(TFseq_minus) <- names(TFmotifs)
for (i in 1:length(TFmotifs)) {
  TFseq_minus[[i]] <- fread(paste(names(TFmotifs[i]), '_hg38seq.bed', sep = ''))
}
#Make sure you still have 400k rows
for (i in 1:length(TFseq_minus)) {
  print(nrow(TFseq_minus[[i]]))
}
#Make all sequences upper case
TFseq_minus = lapply(TFseq_minus, function(x) {toupper(x$V7)})
save(TFseq_minus, file = 'hg38_Minus_motif_sequences.Rdata')
#Break sequences into 5bp chunks (for 5mer)
TFseq_minus = lapply(TFseq_minus, mer.positions, mermask = "NNNNN")
#Find frequency of each possible 5mer for each position
TFseq_minus = lapply(TFseq_minus, position.frequencies, mermask = 'NNNNN')
#Multiply 1/scale_factor by its corresponding kmer frequency for each TF motif. This will produce an
# .rds file for each TF/mask combo (scale factor * k-mer frequency - sfkf)
registerDoParallel(3)
TF_scalefactors_kmerfreq <- vector(mode = "list", length = length(TFseq_minus))
names(TF_scalefactors_kmerfreq) = names(TFseq_minus)

foreach (i = 1:length(TF_scalefactors_kmerfreq)) %dopar% {TF_scalefactors_kmerfreq[[i]] = 
  scalefactor.by.kmerfrequency(scalefactors = scalefactors,
                               kmerfrequency = TFseq_minus[[i]], tfname = names(TFseq_minus[i]))}

#Next we load in each sfkf combo to the same list
TF_sfkf_5mer_minus <- vector(mode = "list", length = length(TFseq_minus))
names(TF_sfkf_5mer_minus) = names(TFseq_minus)


RDSlist = list.files('./')
RDSlist = RDSlist[grep('rm_minus_DNase_', RDSlist)]
RDSlist = RDSlist[grep('.rds', RDSlist)]
for (i in 1:length(RDSlist)) {
  TF_sfkf_5mer_minus[[which(names(TF_sfkf_5mer_minus) == 
          substr(RDSlist[i], 1, nchar(RDSlist[i]) - 48))]][[grep(substr(RDSlist[i],
          nchar(RDSlist[i]) - 40, nchar(RDSlist[i])- 4),
          factorfiles)]] = readRDS(paste( './', RDSlist[i], sep = ''))
          names(TF_sfkf_5mer_minus[[which(names(TF_sfkf_5mer_minus) == 
          substr(RDSlist[i], 1, nchar(RDSlist[i]) - 48))]])[grep(substr(RDSlist[i],
          nchar(RDSlist[i]) - 40, nchar(RDSlist[i])- 4),
          factorfiles)] =  substr(RDSlist[i], 1, nchar(RDSlist[i]) - 4)
          TF_sfkf_5mer_minus[[which(names(TF_sfkf_5mer_minus) == 
          substr(RDSlist[i], 1, nchar(RDSlist[i]) - 48))]][[grep(substr(RDSlist[i],
          nchar(RDSlist[i]) - 40, nchar(RDSlist[i])- 4),
          factorfiles)]] = scaledkmerfreq.sum.minus(TF_sfkf_5mer_minus[[which(names(TF_sfkf_5mer_minus)
          == substr(RDSlist[i], 1,
          nchar(RDSlist[i]) - 48))]][[grep(substr(RDSlist[i], nchar(RDSlist[i]) - 40,
          nchar(RDSlist[i])- 4), factorfiles)]])
}
minus_pre_input = TF_sfkf_5mer_minus

###############################################################################
#Load in unscaled composites
load('DNase_minus_unscaled_compositelist.Rdata')
DNase_unscaled_composites = do.call(rbind, DNase_minus_unscaled_compositelist)
#Trim unscaled composites to the size of predicted area
DNase_unscaled_composites = DNase_unscaled_composites[-c(which(DNase_unscaled_composites$x < -90.5 |
                                                                 DNase_unscaled_composites$x > 90.5)),]
#Normalize unscaled composites to the same value for all composites
DNase_unscaled_composites$est = DNase_unscaled_composites$est/mean(DNase_unscaled_composites$correction)
#Subset out training motifs
DNase_unscaled_composites = DNase_unscaled_composites[which(substr(DNase_unscaled_composites$factor,
                            1,nchar(DNase_unscaled_composites$factor)-9) %in% train),]
#Shift pre predictor variables. Each masks' position shifts by the mask's offset
input_correction = data.frame(matrix(nrow = 32, ncol = 2))
input_correction[1:32,1] = 1:16
input_correction[1:32,2] = 182:197

minus_pred_input = vector(mode = 'list', length = length(minus_pre_input))
names(minus_pred_input) <- names(minus_pre_input)
for (i in 1:length(minus_pre_input)) {
  for (p in 1:length(minus_pre_input[[i]])) {
    minus_pred_input[[i]][[p]] = minus_pre_input[[i]][[p]][input_correction[p,1]:input_correction[p,2]]
  }
  minus_pred_input[[i]] = data.frame(minus_pred_input[[i]])
}

#Assign names
for (i in 1:length(minus_pred_input)) {
  names(minus_pred_input[[i]]) = names(minus_pre_input[[i]])  
}
#Make all mask names the same:
for (i in 1:length(minus_pred_input)) {
  names(minus_pred_input[[i]]) = substr(names(minus_pred_input[[i]]),
            nchar(names(minus_pred_input[[i]]))-36, nchar(names(minus_pred_input[[i]])))
}
#Turn minus_pred_input into DF
minus_pred_input = cbind(seq(-90.5,90.5,1), DNase_unscaled_composites$est, do.call(rbind, minus_pred_input))
colnames(minus_pred_input)[1:2] = c('X', 'unscaled')
minus_pred_input = minus_pred_input[which(minus_pred_input$X >= -39.5 & minus_pred_input$X <= 39.5),]
#This is the plus input for the rule ensemble
save(minus_pred_input, file = 'DNase_TF_minus_pre_input.Rdata')
