library(lattice)
library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
#Load in scale factor files
Benzonase_factor_files = list.files('EnzymeBias_masks')
#Subset only the Benzonase files
Benzonase_factor_files = Benzonase_factor_files[grep('Benzonase', Benzonase_factor_files)]
Benzonase_scale_factors = vector(mode = 'list', length = length(Benzonase_factor_files))
#Determine scale factors for each k-mer position 
for (i in 1:length(Benzonase_factor_files)) {
  Benzonase_scale_factors[[i]] = scalefactor.func(paste('EnzymeBias_masks/',
                                                        Benzonase_factor_files[i], sep = ''))
}
names(Benzonase_scale_factors) = gsub( '_scale_factors.txt', '',Benzonase_factor_files)
names(Benzonase_scale_factors) = gsub( 'Benzonase_', '',names(Benzonase_scale_factors))
#Turn these scale factors into a data.table with an identifier column 
Benzonase_scale_factors_df = rbindlist(Benzonase_scale_factors, idcol = TRUE)
colnames(Benzonase_scale_factors_df)[1] = 'ID'
factor_key = data.frame(levels(factor(Benzonase_scale_factors_df$ID)))
factor_key[,2] = 1:nrow(factor_key)
colnames(factor_key) = c('mask', 'index')
Benzonase_scale_factors_df[, xaxis := factor_key[match(Benzonase_scale_factors_df$ID, factor_key$mask),2]-16.5]
######################################################################################################
#Add benzonase to the list of enzyme scale factors for plotting
Figure1B_list = vector(mode = 'list', length = 1L)
Figure1B_list[[1]] = Benzonase_scale_factors_df[,c(10,14)]
#Combine Benzonase scale factors for plotting
colnames(Figure1B_list[[1]])[1] = 'minusscalefact'
Figure1B_list[[1]] = rbind(Figure1B_list[[1]], Benzonase_scale_factors_df[,c(13,14)])
names(Figure1B_list)[1] = 'Benzonase_scale_factors_df'
######################################################################################################
#Load in scale factor files
Cyanase_factor_files = list.files('EnzymeBias_masks')
#Subset only the Cyanase files
Cyanase_factor_files = Cyanase_factor_files[grep('Cyanase', Cyanase_factor_files)]
Cyanase_scale_factors = vector(mode = 'list', length = length(Cyanase_factor_files))
#Determine scale factors for each k-mer position 
for (i in 1:length(Cyanase_factor_files)) {
  Cyanase_scale_factors[[i]] = scalefactor.func(paste('EnzymeBias_masks/',
                                                      Cyanase_factor_files[i], sep = ''))
}
names(Cyanase_scale_factors) = gsub( '_scale_factors.txt', '',Cyanase_factor_files)
names(Cyanase_scale_factors) = gsub( 'Cyanase_', '',names(Cyanase_scale_factors))
#Turn these scale factors into a data.table with an identifier column 
Cyanase_scale_factors_df = rbindlist(Cyanase_scale_factors, idcol = TRUE)
colnames(Cyanase_scale_factors_df)[1] = 'ID'
factor_key = data.frame(levels(factor(Cyanase_scale_factors_df$ID)))
factor_key[,2] = 1:nrow(factor_key)
colnames(factor_key) = c('mask', 'index')
Cyanase_scale_factors_df[, xaxis := factor_key[match(Cyanase_scale_factors_df$ID, factor_key$mask),2]-16.5]
######################################################################################################
#Add Cyanase to the list of enzyme scale factors for plotting
Figure1B_list[[2]] = Cyanase_scale_factors_df[,c(10,14)]
#Combine Cyanase plus and minus scale factors for plotting:
colnames(Figure1B_list[[2]])[1] = 'minusscalefact'
Figure1B_list[[2]] = rbind(Figure1B_list[[2]], Cyanase_scale_factors_df[,c(13,14)])
names(Figure1B_list)[2] = 'Cyanase_scale_factors_df'
######################################################################################################
######################################################################################################
#Load in scale factor files
MNase_factor_files = list.files('EnzymeBias_masks')
#Subset only the MNase files
MNase_factor_files = MNase_factor_files[grep('MNase', MNase_factor_files)]
MNase_factor_files = MNase_factor_files[grep('scale_factors', MNase_factor_files)]
MNase_scale_factors = vector(mode = 'list', length = length(MNase_factor_files))
#Determine scale factors for each k-mer position 
for (i in 1:length(MNase_factor_files)) {
  MNase_scale_factors[[i]] = scalefactor.func(paste('EnzymeBias_masks/',
                                                    MNase_factor_files[i], sep = ''))
}
names(MNase_scale_factors) = gsub( '_scale_factors.txt', '',MNase_factor_files)
names(MNase_scale_factors) = gsub( 'MNase_', '',names(MNase_scale_factors))
#Turn these scale factors into a data.table with an identifier column 
MNase_scale_factors_df = rbindlist(MNase_scale_factors, idcol = TRUE)
colnames(MNase_scale_factors_df)[1] = 'ID'
factor_key = data.frame(levels(factor(MNase_scale_factors_df$ID)))
factor_key[,2] = 1:nrow(factor_key)
colnames(factor_key) = c('mask', 'index')
MNase_scale_factors_df[, xaxis := factor_key[match(MNase_scale_factors_df$ID, factor_key$mask),2]-16.5]
######################################################################################################
#Add MNase to the list of enzyme scale factors for plotting
Figure1B_list[[3]] = MNase_scale_factors_df[,c(10,14)]
#Combine minus and plus scale factors for plotting
colnames(Figure1B_list[[3]])[1] = 'minusscalefact'
Figure1B_list[[3]] = rbind(Figure1B_list[[3]], MNase_scale_factors_df[,c(13,14)])
names(Figure1B_list)[3] = 'MNase_scale_factors_df'
######################################################################################################
######################################################################################################
#Load in scale factor files
DNase_factor_files = list.files('EnzymeBias_masks')
#Subset only the DNase files
DNase_factor_files = DNase_factor_files[grep('DNase', DNase_factor_files)]
DNase_scale_factors = vector(mode = 'list', length = length(DNase_factor_files))
#Determine scale factors for each k-mer position 
for (i in 1:length(DNase_factor_files)) {
  DNase_scale_factors[[i]] = scalefactor.func(paste('EnzymeBias_masks/',
                                                    DNase_factor_files[i], sep = ''))
}
names(DNase_scale_factors) = gsub( '_scale_factors.txt', '',DNase_factor_files)
names(DNase_scale_factors) = gsub( 'DNase_', '',names(DNase_scale_factors))
#Turn these scale factors into a data.table with an identifier column 
DNase_scale_factors_df = rbindlist(DNase_scale_factors, idcol = TRUE)
colnames(DNase_scale_factors_df)[1] = 'ID'
factor_key = data.frame(levels(factor(DNase_scale_factors_df$ID)))
factor_key[,2] = 1:nrow(factor_key)
colnames(factor_key) = c('mask', 'index')
DNase_scale_factors_df[, xaxis := factor_key[match(DNase_scale_factors_df$ID, factor_key$mask),2]-16.5]
######################################################################################################
#Add DNase to the list of enzyme scale factors for plotting
Figure1B_list[[4]] = DNase_scale_factors_df[,c(10,14)]
#Combine minus and plus scale factors for plotting
colnames(Figure1B_list[[4]])[1] = 'minusscalefact'
Figure1B_list[[4]] = rbind(Figure1B_list[[4]], DNase_scale_factors_df[,c(13,14)])
names(Figure1B_list)[4] = 'DNase_scale_factors_df'
######################################################################################################
######################################################################################################
#Load in scale factor files
Tn5_factor_files = list.files('EnzymeBias_masks')
#Subset only the Tn5 files with masks we desire
Tn5_factor_files = Tn5_factor_files[grep('Tn5', Tn5_factor_files)]
Tn5_factor_files = Tn5_factor_files[(grep('XXXXXXXXXXNNNNNXXXXXXXXCXXXXXXXXXXXXXXXXXXXXXXX', Tn5_factor_files)):
                                      (grep('XXXXXXXXXXNNNNNXXXXXXXXCXXXXXXXXXXXXXXXXXXXXXXX', Tn5_factor_files)+30)]
Tn5_scale_factors = vector(mode = 'list', length = length(Tn5_factor_files))
#Determine scale factors for each k-mer position 
for (i in 1:length(Tn5_factor_files)) {
  Tn5_scale_factors[[i]] = scalefactor.func(paste('EnzymeBias_masks/',
                                                  Tn5_factor_files[i], sep = ''))
}
names(Tn5_scale_factors) = gsub( '_scale_factors.txt', '',Tn5_factor_files)
names(Tn5_scale_factors) = gsub( 'Tn5_', '',names(Tn5_scale_factors))
#Turn these scale factors into a data.table with an identifier column 
Tn5_scale_factors_df = rbindlist(Tn5_scale_factors, idcol = TRUE)
colnames(Tn5_scale_factors_df)[1] = 'ID'
factor_key = data.frame(levels(factor(Tn5_scale_factors_df$ID)))
factor_key[,2] = 1:nrow(factor_key)
colnames(factor_key) = c('mask', 'index')
Tn5_scale_factors_df[, xaxis := factor_key[match(Tn5_scale_factors_df$ID, factor_key$mask),2]-16]
######################################################################################################
#Add Tn5 to the list of enzyme scale factors for plotting
Figure1B_list[[5]] = Tn5_scale_factors_df[,c(10,14)]
#Combine minus and plus scale factors for plotting
colnames(Figure1B_list[[5]])[1] = 'minusscalefact'
Figure1B_list[[5]] = rbind(Figure1B_list[[5]], Tn5_scale_factors_df[,c(13,14)])
names(Figure1B_list)[5] = 'Tn5_scale_factors_df'
######################################################################################################
save(Figure1B_list, file = 'Figure1B_scalefactor_list.Rdata')
