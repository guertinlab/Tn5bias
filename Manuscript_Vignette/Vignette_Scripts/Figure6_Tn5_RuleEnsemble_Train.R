library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
library('pre')
library('ggplot2')
options(scipen = 100)
################################################################################
train = c( 'ASCL1', 'ATF3', 'CLOCK', 'CTCF',
           'DLX2', 'EBF1', 'EGR1', 'ELF1',
           'FOS', 'FOXA1', 'FOXK2', 'GATA2',
           'GLIS1', 'IRF1', 'JUN', 'LEF1',
           'MEIS2', 'MLX', 'MYC', 'PPARG', 
           'RUNX1', 'SPIB', 'SRY', 'STAT1', 'TGIF1')
test = c('AR', 'CEBPB', 'E2F1', 'ESR1',
         'FERD3L', 'HOXC12', 'HSF1', 'MAX',
         'MEF2A', 'MGA', 'NFATC3', 'NR2F2',
         'POU3F1', 'REST', 'Six3', 'SP1',
         'USF1', 'TEAD1')

#Load in 5mer Rules Ensemble input and add other inputs from here:
load('Tn5_TF_Plus_pre_input_5mers.Rdata')
Plus_pred_input = Plus_pred_input[order(names(Plus_pred_input))]
#Remove test data
Plus_pred_input = Plus_pred_input[which(substr(names(Plus_pred_input),
                                               1, nchar(names(Plus_pred_input))-8) %in% train)]
#Trim 5mer RE input to match smaller widths:
Plus_pre_input = vector(mode = 'list', length = length(Plus_pred_input))
for (i in 1:length(Plus_pred_input)) {
  for (p in 1:length(Plus_pred_input[[i]])) {
    Plus_pre_input[[i]][[p]] = Plus_pred_input[[i]][[p]][c(which(Plus_pred_input[[i]][[p]]$x >= -73 &
                                                                   Plus_pred_input[[i]][[p]]$x <= 57)),1]
    names(Plus_pre_input)[i] = names(Plus_pred_input[i])
    Plus_pre_input[[i]][[p]] = as.data.frame(Plus_pre_input[[i]][[p]])
    
  }
  Plus_pre_input[[i]] = as.data.frame(Plus_pre_input[[i]])
  names(Plus_pre_input[[i]]) = names(Plus_pred_input[[i]])
}
#load in unscaled composites and reduce the size to the predicted 
#area 
load('../data/Tn5_plus_unscaled_RM_composites.Rdata')
Tn5_unscaled_composites = 
  Tn5_unscaled_composites[which(substr(Tn5_unscaled_composites$factor,
                                       1, nchar(Tn5_unscaled_composites$factor)-8) %in% train),]
Tn5_unscaled_composites = 
  Tn5_unscaled_composites[-c(which(Tn5_unscaled_composites$x < - 73 |  Tn5_unscaled_composites$x > 57)),]
Tn5_pre_data_plus = cbind(Tn5_unscaled_composites, do.call(rbind, Plus_pre_input))
rm(Plus_pre_input, Plus_pred_input)
#Next, load in 6mers...
load('../RE_input/Tn5_TF_Plus_pre_input_6mers.Rdata')
Plus_pred_input = Plus_pred_input[order(names(Plus_pred_input))]
Plus_pred_input = 
  Plus_pred_input[which(substr(names(Plus_pred_input), 1, nchar(names(Plus_pred_input))-8) %in% train)]
#Trim 6mer RE input to match smaller widths:
Plus_pre_input = vector(mode = 'list', length = length(Plus_pred_input))
for (i in 1:length(Plus_pred_input)) {
  for (p in 1:length(Plus_pred_input[[i]])) {
    Plus_pre_input[[i]][[p]] = Plus_pred_input[[i]][[p]][c(which(Plus_pred_input[[i]][[p]]$x >= -73 &
                                                                   Plus_pred_input[[i]][[p]]$x <= 57)),1]
    names(Plus_pre_input)[i] = names(Plus_pred_input[i])
    Plus_pre_input[[i]][[p]] = as.data.frame(Plus_pre_input[[i]][[p]])
  }
  Plus_pre_input[[i]] = as.data.frame(Plus_pre_input[[i]])
  names(Plus_pre_input[[i]]) = names(Plus_pred_input[[i]])
}
#Append 6mers:
Tn5_pre_data_plus = cbind(Tn5_pre_data_plus, do.call(rbind, Plus_pre_input))
rm(Plus_pre_input, Plus_pred_input)
#Next, load in 7mers...
load('../RE_input/Tn5_TF_Plus_pre_input_7mers.Rdata')
Plus_pred_input = Plus_pred_input[order(names(Plus_pred_input))]
Plus_pred_input = Plus_pred_input[which(substr(names(Plus_pred_input),
                                               1, nchar(names(Plus_pred_input))-8) %in% train)]
#Trim 7mer RE input to match smaller widths:
Plus_pre_input = vector(mode = 'list', length = length(Plus_pred_input))
for (i in 1:length(Plus_pred_input)) {
  for (p in 1:length(Plus_pred_input[[i]])) {
    Plus_pre_input[[i]][[p]] = Plus_pred_input[[i]][[p]][c(which(Plus_pred_input[[i]][[p]]$x >= -73 &
                                                                   Plus_pred_input[[i]][[p]]$x <= 57)),1]
    names(Plus_pre_input)[i] = names(Plus_pred_input[i])
    Plus_pre_input[[i]][[p]] = as.data.frame(Plus_pre_input[[i]][[p]])
  }
  Plus_pre_input[[i]] = as.data.frame(Plus_pre_input[[i]])
  names(Plus_pre_input[[i]]) = names(Plus_pred_input[[i]])
}
#Append 7mers:
Tn5_pre_data_plus = cbind(Tn5_pre_data_plus, do.call(rbind, Plus_pre_input))
rm(Plus_pre_input, Plus_pred_input)
#Next, load in spaced 3mers...
load('../RE_input/Tn5_TF_Plus_pre_input_spaced3mers.Rdata')
Plus_pred_input = Plus_pred_input[order(names(Plus_pred_input))]
Plus_pred_input = Plus_pred_input[which(substr(names(Plus_pred_input),
                                               1, nchar(names(Plus_pred_input))-8) %in% train)]
#Trim 7mer RE input to match smaller widths:
Plus_pre_input = vector(mode = 'list', length = length(Plus_pred_input))
for (i in 1:length(Plus_pred_input)) {
  for (p in 1:length(Plus_pred_input[[i]])) {
    Plus_pre_input[[i]][[p]] = Plus_pred_input[[i]][[p]][c(which(Plus_pred_input[[i]][[p]]$x >= -73 &
                                                                   Plus_pred_input[[i]][[p]]$x <= 57)),1]
    names(Plus_pre_input)[i] = names(Plus_pred_input[i])
    Plus_pre_input[[i]][[p]] = as.data.frame(Plus_pre_input[[i]][[p]])
    
  }
  Plus_pre_input[[i]] = as.data.frame(Plus_pre_input[[i]])
  names(Plus_pre_input[[i]]) = names(Plus_pred_input[[i]])
}
#Append spaced 3mers:
Tn5_pre_data_plus = cbind(Tn5_pre_data_plus, do.call(rbind, Plus_pre_input))
rm(Plus_pre_input, Plus_pred_input, Tn5_unscaled_composites)
##Repeat this process for minus strand data
#Load in original 5mer Rules Ensemble input and add other inputs from here:
load('../RE_input/Tn5_TF_minus_pre_input_5mers.Rdata')
minus_pred_input = minus_pred_input[order(names(minus_pred_input))]
#Remove test data
minus_pred_input = minus_pred_input[which(substr(names(minus_pred_input),
                                                 1, nchar(names(minus_pred_input))-9) %in% train)]
#Trim 5mer RE input to match smaller widths:
minus_pre_input = vector(mode = 'list', length = length(minus_pred_input))
for (i in 1:length(minus_pred_input)) {
  for (p in 1:length(minus_pred_input[[i]])) {
    minus_pre_input[[i]][[p]] = minus_pred_input[[i]][[p]][c(which(minus_pred_input[[i]][[p]]$x >= -73 &
                                                                     minus_pred_input[[i]][[p]]$x <= 57)),1]
    names(minus_pre_input)[i] = names(minus_pred_input[i])
    minus_pre_input[[i]][[p]] = as.data.frame(minus_pre_input[[i]][[p]])
  }
  minus_pre_input[[i]] = as.data.frame(minus_pre_input[[i]])
  names(minus_pre_input[[i]]) = names(minus_pred_input[[i]])
}

#load in unscaled composites and reduce the size to the predicted 
#area (last 4 bases are unmappable because of 5mer)
load('../data/Tn5minus_unscaled_compositelist.Rdata')
Tn5minus_unscaled_composites = do.call(rbind, Tn5minus_unscaled_compositelist)
Tn5minus_unscaled_composites$est = Tn5minus_unscaled_composites$est/0.00309736

Tn5minus_unscaled_composites = Tn5minus_unscaled_composites[which(substr(Tn5minus_unscaled_composites$factor,
                                          1, nchar(Tn5minus_unscaled_composites$factor)-9) %in% train),]
Tn5minus_unscaled_composites = Tn5minus_unscaled_composites[-c(which(Tn5minus_unscaled_composites$x < - 73 |
                                                           Tn5minus_unscaled_composites$x > 57)),]

Tn5_pre_data_minus = cbind(Tn5minus_unscaled_composites, do.call(rbind, minus_pre_input))
rm(minus_pre_input, minus_pred_input, Tn5minus_unscaled_compositelist, Tn5minus_unscaled_composites)
#Next, load in 6mers...
load('../RE_input/Tn5_TF_minus_pre_input_6mers.Rdata')
minus_pred_input = minus_pred_input[order(names(minus_pred_input))]
minus_pred_input = minus_pred_input[which(substr(names(minus_pred_input),
                                                 1, nchar(names(minus_pred_input))-9) %in% train)]
#Trim 6mer RE input to match smaller widths:
minus_pre_input = vector(mode = 'list', length = length(minus_pred_input))
for (i in 1:length(minus_pred_input)) {
  for (p in 1:length(minus_pred_input[[i]])) {
    minus_pre_input[[i]][[p]] = minus_pred_input[[i]][[p]][c(which(minus_pred_input[[i]][[p]]$x >= -73 &
                                                                     minus_pred_input[[i]][[p]]$x <= 57)),1]
    names(minus_pre_input)[i] = names(minus_pred_input[i])
    minus_pre_input[[i]][[p]] = as.data.frame(minus_pre_input[[i]][[p]])
    
  }
  minus_pre_input[[i]] = as.data.frame(minus_pre_input[[i]])
  names(minus_pre_input[[i]]) = names(minus_pred_input[[i]])
}
#Append 6mers:
Tn5_pre_data_minus = cbind(Tn5_pre_data_minus, do.call(rbind, minus_pre_input))
rm(minus_pre_input, minus_pred_input)
#Next, load in 7mers...
load('../RE_input/Tn5_TF_minus_pre_input_7mers.Rdata')
minus_pred_input = minus_pred_input[order(names(minus_pred_input))]
minus_pred_input = minus_pred_input[which(substr(names(minus_pred_input), 1,
                                                 nchar(names(minus_pred_input))-9) %in% train)]
#Trim 7mer RE input to match smaller widths:
minus_pre_input = vector(mode = 'list', length = length(minus_pred_input))
for (i in 1:length(minus_pred_input)) {
  for (p in 1:length(minus_pred_input[[i]])) {
    minus_pre_input[[i]][[p]] = minus_pred_input[[i]][[p]][c(which(minus_pred_input[[i]][[p]]$x >= -73 &
                                                                     minus_pred_input[[i]][[p]]$x <= 57)),1]
    names(minus_pre_input)[i] = names(minus_pred_input[i])
    minus_pre_input[[i]][[p]] = as.data.frame(minus_pre_input[[i]][[p]])
    
  }
  minus_pre_input[[i]] = as.data.frame(minus_pre_input[[i]])
  names(minus_pre_input[[i]]) = names(minus_pred_input[[i]])
}
#Append 7mers:
Tn5_pre_data_minus = cbind(Tn5_pre_data_minus, do.call(rbind, minus_pre_input))
rm(minus_pre_input, minus_pred_input)
#Next, load in spaced 3mers...
load('../RE_input/Tn5_TF_minus_pre_input_spaced3mers.Rdata')
minus_pred_input = minus_pred_input[order(names(minus_pred_input))]
minus_pred_input = minus_pred_input[which(substr(names(minus_pred_input),
                                                 1, nchar(names(minus_pred_input))-9) %in% train)]
#Trim spaced 3mer RE input to match smaller widths:
minus_pre_input = vector(mode = 'list', length = length(minus_pred_input))
for (i in 1:length(minus_pred_input)) {
  for (p in 1:length(minus_pred_input[[i]])) {
    minus_pre_input[[i]][[p]] = minus_pred_input[[i]][[p]][c(which(minus_pred_input[[i]][[p]]$x >= -73 &
                                                                     minus_pred_input[[i]][[p]]$x <= 57)),1]
    names(minus_pre_input)[i] = names(minus_pred_input[i])
    minus_pre_input[[i]][[p]] = as.data.frame(minus_pre_input[[i]][[p]])
    
  }
  minus_pre_input[[i]] = as.data.frame(minus_pre_input[[i]])
  names(minus_pre_input[[i]]) = names(minus_pred_input[[i]])
}
#Append spaced 3mers:
Tn5_pre_data_minus = cbind(Tn5_pre_data_minus, do.call(rbind, minus_pre_input))
rm(minus_pre_input, minus_pred_input)
###Combine minus and plus data together:
Tn5_RE_input = rbind(Tn5_pre_data_plus, Tn5_pre_data_minus)
##Trim off beginning 7 positions:
RE_colnames = strsplit(colnames(Tn5_RE_input)[6:948], '')
RE_colnames = as.data.frame(do.call(rbind, RE_colnames))
REmoval_columns = which(RE_colnames[,1] == 'N' | RE_colnames[,2] == 'N' |
                          RE_colnames[,3] == 'N' | RE_colnames[,4] == 'N'|
                          RE_colnames[,5] == 'N' | RE_colnames[,6] == 'N' |
                          RE_colnames[,7] == 'N') + 5
Tn5_RE_input = as.data.frame(Tn5_RE_input)
Tn5_RE_input = Tn5_RE_input[,-c(REmoval_columns)]
Tn5_RE_input = as.data.table(Tn5_RE_input)



#Remove columns 1,3,4,5
Tn5_RE_input = Tn5_RE_input[,-c(1,3,4,5)]
sapply(Tn5_RE_input, function(x) sum(is.na(x)))
################################################################################
#Make linear model
set.seed(42)
Tn5_linear_model = pre(est ~  ., data = Tn5_RE_input, type = 'linear', sampfrac = 0.5,
                       verbose = TRUE, ntrees = 50000,
                       use.grad = FALSE, nfolds = 5, winsfrac = 0, maxdepth = 3L,
                       learnrate = 0.1, intercept = FALSE)
save(Tn5_linear_model, file = 'Tn5_Linear_model.Rdata')
load('Tn5_Linear_model.Rdata')
################################################################################
#print(Tn5_pre_model, penalty.par.val = Tn5_pre_model$glmnet.fit$lambda.min)
#Plot variable importances for PRE model:
imps = importance(Tn5_linear_model, penalty.par.val =
                    Tn5_linear_model$glmnet.fit$lambda.min, abbreviate = FALSE)
base_imps = imps$baseimps
##Take top 10% of masks (94) for Rules Ensemble:
RE_masks = base_imps$rule[1:66]
write.table(RE_masks, file = 'RulesEnsemble_top10_masks.txt', row.names = FALSE,
            quote = FALSE, sep = '\n', col.names = FALSE)

#Rm all of importance dfs
rm(imps, base_imps)
###Create data frame with only top 10%:
Tn5_RE_input = as.data.frame(Tn5_RE_input)
Tn5_RE_train = Tn5_RE_input[,c(1, which(colnames(Tn5_RE_input) %in% RE_masks))]
save(Tn5_RE_train, file = 'Tn5_RE_input_top10masks.Rdata')
#Make Prediction using linear+rules PRE model
set.seed(42)
Tn5_RE_model = pre(est ~  ., data = Tn5_RE_train, type = 'both',
                   sampfrac = 0.5, verbose = TRUE, ntrees = 100000,
                   use.grad = FALSE, nfolds = 10, winsfrac = 0,
                   maxdepth = 3L, learnrate = 0.1, intercept = FALSE)

save(Tn5_RE_model, file = 'Tn5_RE_model.Rdata')
load('Tn5_RE_model.Rdata')

#Write tables for each desired model coefficients/rules
print(Tn5_RE_model, penalty.par.val = Tn5_RE_model$glmnet.fit$lambda.min)
Tn5_RE_model_coef <- coef(Tn5_RE_model, penalty.par.val = Tn5_RE_model$glmnet.fit$lambda.min)
Tn5_RE_model_coef = Tn5_RE_model_coef[which(Tn5_RE_model_coef$coefficient != 0),]
write.table(Tn5_RE_model_coef[,2:3], file = 'Tn5_Rules_Ensemble_model.txt',
            quote = FALSE, sep = ',', row.names = FALSE)

#Plot Importances
imps = importance(Tn5_RE_model, penalty.par.val = Tn5_RE_model$glmnet.fit$lambda.min, abbreviate = FALSE)
var_imps = imps$varimps
var_imps$varname = as.factor(var_imps$varname)
redundant_var_imps = redundantpos_nml(var_imps, labelnum = 1, impnum = 2)
redundant_var_imps = redundant_var_imps[-c(1:7),]
redundant_var_imps = as.data.frame(redundant_var_imps)

##Shift positions by 4
imppos = strsplit(redundant_var_imps$position, '')
imppos = as.data.frame(do.call(rbind, imppos))
Ccol = imppos[,which(imppos[1,] == 'C')]
Ccolnum = which(imppos[1,] == 'C')

imppos[,48:51] = 'X'
imppos = imppos[,c(48:51,1:47)]
colnames(imppos) = paste('V',1:ncol(imppos), sep = '')
Ccol = imppos[,which(imppos[1,] == 'C')]
Ccolnum = which(imppos[1,] == 'C')
imppos = imppos[,-c(which(imppos[1,] == 'C'))]

for (i in 1:nrow(imppos)) {
  Npos = which(imppos[i,] == 'N')
  imppos[i,Npos] = 'X'
  imppos[i,(Npos-4)] = 'N'
}

imppos = cbind(imppos, Ccol)
imppos = imppos[,c(1:(Ccolnum-1),ncol(imppos),(Ccolnum):(ncol(imppos)-1))]

for (i in 1:nrow(imppos)) {
  imppos[i,] = paste(imppos[i,], collapse = '')
}
imppos = imppos[,1]

redundant_var_imps[,2] = imppos
redundant_var_imps$position = seq(-19,19,1) 
save(redundant_var_imps, file = 'Tn5_RE_redundant_var_imps.Rdata')

#Plot redundant_var_imps
pdf(file = 'Figure6A_Redundant_Tn5_RE_variable_importances.pdf', width = 12, height = 12)
impbarplot <-ggplot(data=redundant_var_imps, aes(x=position, y=possum)) 
impbarplot +
  labs(x= 'Position Relative to Central Base', y= 'Importance') +
  geom_bar(stat="identity", fill="black") +
  theme_classic() +
  theme(axis.title.x = element_text(size=52, colour = "black", face = 'bold'),
        axis.title.y = element_text(size=52, colour = "black", face = 'bold'),
        axis.text.x = element_text(size=38, colour = "black"),
        axis.text.y = element_text(size=38, colour = "black"),
        axis.line = element_line(colour = 'black', size = 1),
        plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(vjust=1), plot.margin = margin(1,1,1.5,1.2, "cm"),
        axis.ticks.length=unit(.75, "cm"), axis.ticks=element_line(size=2)) +
  scale_x_continuous(breaks = c(-15,0,15))
dev.off()
