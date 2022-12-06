library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
################################################################################
#Load in unscaled bedgraph in order to get read depth
system('bigWigtoBedGraph ../Figure1/C1_gDNA_rep1.bigWig C1_gDNA_rep1_unscaled.bedGraph')
unscaled_bed = fread('C1_gDNA_rep1_unscaled.bedGraph')
unscaled_read_depth = sum(unscaled_bed$V4)
rm(unscaled_bed)
#Load in bedgraph which contains all input mask bedgraphs
x <- fread('RE_masks/Tn5_top10_union.bedGraph')
#Load in column names for the above bedgraph
masknames = fread('RE_masks/Tn5_top10_names.txt', header = FALSE)
masknames = masknames$V1
masknames = substr(masknames, 5, nchar(masknames))
masknames = c('chr', 'start', 'stop', masknames)
colnames(x) = masknames
#Import rule ensemble model

#Import the formula from the rule ensemble 
prerules = import.rules('Tn5_Rules_Ensemble_model.txt')
prerules = gsub(' + ', ' + \n', prerules, fixed = TRUE)
#Format the formula into a function and write a file containing this function
write.table('', file = 'Tn5_RuleEnsemble_model.R')
fileConn<-file("Tn5_RuleEnsemble_model.R")
writeLines(c('RE_scale_Tn5 = function (input) { input[, RuleEnsemble :=', prerules, '] \nreturn(input)}'), fileConn)
close(fileConn)
source('Tn5_RuleEnsemble_model.R')
#This function can also be called from GitHub:
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_RuleEnsemble_model.R')
#Apply the model to Tn5 bedGraph
x = RE_scale_Tn5(x)
#Get model read depth
pre_read_depth = sum(x$RuleEnsemble)
#Create a ratio to scale Rule ensemble read depth to unscaled read depth
RDS = unscaled_read_depth/pre_read_depth
x[, pre_RDS := RuleEnsemble*RDS ]
#Write the scaled DNase data in bedGraph format
write.table(x[,c(1:3,71)], file = 'Tn5_Rule_Ensemble_Scaled.bedGraph', col.names = FALSE, 
            row.names=FALSE, sep = '\t', quote=FALSE)
#Convert the bedGraph to bigWig format
system('bedGraphToBigWig Tn5_Rule_Ensemble_Scaled.bedGraph ../Figure1/hg38.fa.chrom.sizes Tn5_Rule_Ensemble_Scaled.bigWig')