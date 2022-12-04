library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)

#Import kmer counts and compute scale factors for each....
factorfiles = list.files('Tn5_kmercount/')
factorfiles = factorfiles[grep('_scale_factors.txt', factorfiles)]
gsub_factorfiles = gsub('C', '', factorfiles)
factorfiles = factorfiles[-grep('NNNNN', gsub_factorfiles)]

unique_factorfiles_spacing = substr(factorfiles, 5, nchar(factorfiles) -18)
unique_factorfiles_spacing = gsub('C', '', unique_factorfiles_spacing)
unique_factorfiles_spacing = strsplit(unique_factorfiles_spacing, 'NNN')
unique_factorfiles = NULL
for (i in 1:length(unique_factorfiles_spacing)) {
  unique_factorfiles[i] = unique_factorfiles_spacing[[i]][2]
}
unique_factorfiles_spacing = gsub('^X', 'NNNX', unique_factorfiles)
unique_factorfiles_spacing = gsub('X$', 'XNNN', unique_factorfiles_spacing)
unique_factorfiles_spacing = unique(unique_factorfiles_spacing)
rm(unique_factorfiles,gsub_factorfiles)

scalefactors <- vector(mode = "list", length = length(factorfiles))
names(scalefactors) <- factorfiles

for (i in 1:length(factorfiles)) {
  scalefactors[[i]] <- scalefactor.func(paste('Tn5_kmercount/',
                                              factorfiles[i], sep = ''))
}
#Want scalefactor masks to have Xs instead of Ns
for (i in 1:length(scalefactors)) {
  scalefactors[[i]]$V2  = gsub('N', 'X', scalefactors[[i]]$V2)
}
#Match scalefactor with its kmerfreq
scalefactor_names = substr(names(scalefactors), 5, 51)
scalefactor_names = gsub('C', '', scalefactor_names)

TFseq_files = list.files('spaced_3mer/plus')
#Here we make a 'key' for each spacing of 3mers. 
#Each key lists the spacings of k-mer frequency to be
#multiplied by their respective scale factors
system('mkdir HPC_scripts')
setwd('HPC_scripts')
for (i in 1:length(TFseq_files)) {
  sfkf_key = names(scalefactors)[grep(substr(TFseq_files[i], 1,
                                             nchar(TFseq_files[i])-26), scalefactor_names)]
  save(sfkf_key, file = paste(substr(TFseq_files[i], 1,
                                     nchar(TFseq_files[i])-26), '_sfkf_key.Rdata', sep = ''))
}
#Save all scalefactors
save(scalefactors, file = 'Tn5_spaced6mer_scalefactors.Rdata')
#This section writes a script for each possible spacing of 3-mers within the window
system('mkdir HPC_scripts')
setwd('HPC_scripts')
for (i in 1:length(TFseq_files)) {
  writeLines(c("library(data.table)",
               "source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')", 
               "library('parallel')",
               "library('foreach')", 
               "library('doParallel')",
               paste("load('spaced_6mer/plus/", TFseq_files[i], "')", sep = ''),
               paste("load('", paste(substr(TFseq_files[i], 
                                                         1, nchar(TFseq_files[i])-26),
                                                  '_sfkf_key.Rdata', sep = ''), "')", sep = ''),
               "load('Tn5_spaced6mer_scalefactors.Rdata')",
               "setwd('../plus_Rivanna_output')",
               "options(scipen = 100)",
               "registerDoParallel(20)",
               "TF_scalefactors_kmerfreq <- vector(mode = 'list', length = length(TFseq_plus_s))",
               "names(TF_scalefactors_kmerfreq) = names(TFseq_plus_s)",
               "foreach (i = 1:length(TF_scalefactors_kmerfreq)) %dopar% {TF_scalefactors_kmerfreq[[i]] = scalefactor.by.kmerfrequency(scalefactors = scalefactors[which(names(scalefactors) %in% sfkf_key)],
             kmerfrequency = TFseq_plus_s[[i]], tfname = names(TFseq_plus_s[i]))}",
               "rm(TFseq_plus_s)",
               "setwd('../minus_Rivanna_output')",
               paste("load('spaced_6mer/minus/", gsub('plus_', 'minus_', TFseq_files[i]), "')", sep = ''),
               "TF_scalefactors_kmerfreq <- vector(mode = 'list', length = length(TFseq_minus_s))",
               "names(TF_scalefactors_kmerfreq) = names(TFseq_minus_s)",
               "foreach (i = 1:length(TF_scalefactors_kmerfreq)) %dopar% {TF_scalefactors_kmerfreq[[i]] = scalefactor.by.kmerfrequency(scalefactors = scalefactors[which(names(scalefactors) %in% sfkf_key)],
             kmerfrequency = TFseq_minus_s[[i]], tfname = paste(names(TFseq_minus_s[i]), sep = ''))}"
  ),
  paste(substr(TFseq_files[i], 1, nchar(TFseq_files[i])-26), "_sfkf_Rivanna.R", sep = ''))
}
#This section writes a slurm script to run each of the above scripts using an HPC
for (i in 1:length(TFseq_files)) {
  writeLines(c("#!/bin/bash",
               "#SBATCH --time=3:00:00",
               "#SBATCH -A yourlab",
               "#SBATCH -n 20",
               "#SBATCH -p standard",
               "#SBATCH --mem=9000",
               paste("#SBATCH -o Tn5_", substr(TFseq_files[i],
                                               1, nchar(TFseq_files[i])-26), "_sfkf.out", sep = ''),
               "module load gcc/9.2.0",
               "module load intel/18.0",
               "module load intelmpi/18.0",
               "module load R/4.0.3",
               paste("Rscript ", substr(TFseq_files[i], 1,
                                        nchar(TFseq_files[i])-26), "_sfkf_Rivanna.R", sep = "")),
             paste(substr(TFseq_files[i], 1,
                          nchar(TFseq_files[i])-26), "_sfkf_Rivanna.slurm", sep = ''))
}