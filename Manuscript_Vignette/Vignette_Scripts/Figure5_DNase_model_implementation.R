library(data.table)
options(scipen = 100)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
################################################################################
#Import union bedGraph and set column names
x <- fread('DNase_5mer_union.bedGraph')
masknames = fread('DNase_allmasks_names.txt', header = FALSE)
masknames = masknames$V1
masknames = substr(masknames, 7, nchar(masknames))
masknames = c('chr', 'start', 'stop', masknames)
colnames(x) = masknames
#Get the unscaled read depth for later scaling
unscaled_read_depth = sum(x$Naked_unscaled)
#Import the formula from the rule ensemble 
prerules = import.rules('DNase_pre_coeff.txt')
prerules = gsub(' + ', ' + \n', prerules, fixed = TRUE)
#Format the formula into a function and write a file containing this function
write.table('', file = 'DNase_RuleEnsemble_model.R')
fileConn<-file("DNase_RuleEnsemble_model.R")
writeLines(c('RE_scale_DNase = function (input) { input[, RuleEnsemble :=', prerules, '] \nreturn(input)}'), fileConn)
close(fileConn)
source('DNase_RuleEnsemble_model.R')
#This function can also be called from GitHub:
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/DNase_RuleEnsemble_model.R')
#Apply the model to the DNase data
RE_scale_DNase(x)
#Get model read depth
pre_read_depth = sum(x$RuleEnsemble)
#Create a ratio to scale Rule ensemble read depth to unscaled read depth
RDS = unscaled_read_depth/pre_read_depth
x[, pre_RDS := RuleEnsemble*RDS ]
#Write the scaled DNase data in bedGraph format
write.table(x[,c(1:3,24)], file = 'DNase_Rule_Ensemble_Scaled.bedGraph', col.names = FALSE, 
            row.names=FALSE, sep = '\t', quote=FALSE)
#Convert the bedGraph to bigWig format
system('bedGraphToBigWig DNase_Rule_Ensemble_Scaled.bedGraph ../Figure1/hg38.fa.chrom.sizes DNase_Rule_Ensemble_Scaled.bigWig')
