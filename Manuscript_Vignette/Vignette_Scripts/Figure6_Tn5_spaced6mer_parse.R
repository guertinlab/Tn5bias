library(data.table)
options(scipen = 100)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')

#Load plus spaced 6-mer input 
load('Tn5_TF_plus_pre_input_spaced6mers.Rdata')
#Trim spaced 6-mer input to the lowest row number:
kmerfreq_rows = NULL
for (i in 1:length(TF_plus_scalefactors_kmerfreq)) {
  kmerfreq_rows[i] = nrow(TF_plus_scalefactors_kmerfreq[[1]][[i]])
}
min(kmerfreq_rows)

Plus_pre_input = TF_plus_scalefactors_kmerfreq
for (i in 1:length(Plus_pre_input)) {
  for (p in 1:length(Plus_pre_input[[i]])) {
    Plus_pre_input[[i]][[p]] = Plus_pre_input[[i]][[p]][1:min(kmerfreq_rows)]
  }
}
#Shift each mask by the farthest upstream 3-mer's offset:
itnum = NULL
for (i in 1:40) {
  itnum[i] = 41-i
}

input_correction = data.frame(matrix(ncol = 2))
correction_store = data.frame(matrix(ncol = 2))
for (i in 1:length(itnum)) {
  correction_store[1:itnum[i],1] = rep(1+(i-1), itnum[i]) 
  correction_store[1:itnum[i],2] = rep(137+(i-1), itnum[i]) 
  input_correction = rbind(input_correction, correction_store)
  correction_store = data.frame(matrix(ncol = 2))
}
input_correction = input_correction[-c(1),]
rownames(input_correction) = 1:nrow(input_correction)
#Each mask shifts predictive output by the shift relative to the cut site
#Correct this by assigning the values of each window for shifting
Plus_pred_input = vector(mode = 'list', length = length(Plus_pre_input))
names(Plus_pred_input) <- names(Plus_pre_input)
for (i in 1:length(Plus_pre_input)) {
  for (p in 1:length(Plus_pre_input[[i]])) {
    Plus_pred_input[[i]][[p]] = 
      Plus_pre_input[[i]][[p]][input_correction[p,1]:input_correction[p,2]]
  }
  Plus_pred_input[[i]] = data.frame(Plus_pred_input[[i]])
}
#Assign names
for (i in 1:length(Plus_pred_input)) {
  names(Plus_pred_input[[i]]) = names(Plus_pre_input[[i]])  
}

for (i in 1:length(Plus_pred_input)) {
  names(Plus_pred_input[[i]]) = substr(names(Plus_pred_input[[i]]),
       nchar(names(Plus_pred_input[[i]])) - 54, nchar(names(Plus_pred_input[[i]]))-8)
}

###Add x axis to input for future alignment!
for (i in 1:length(Plus_pred_input)) {
  for (p in 1:length(Plus_pred_input[[i]])) {
    Plus_pred_input[[i]][[p]] = as.data.frame(Plus_pred_input[[i]][[p]])
    Plus_pred_input[[i]][[p]][,2] = seq(-73,min(kmerfreq_rows)-73,1)
    colnames(Plus_pred_input[[i]][[p]]) = c('sfkf', 'x')
  }
}

save(Plus_pred_input, file = 'Tn5_TF_plus_pre_input_spaced6mers.Rdata')

#Load Minus spaced 6-mer input 
load('Tn5_TF_Minus_pre_input_spaced6mers.Rdata')
#Trim spaced 6-mer input to the lowest row number:
kmerfreq_rows = NULL
for (i in 1:length(TF_Minus_scalefactors_kmerfreq)) {
  kmerfreq_rows[i] = nrow(TF_Minus_scalefactors_kmerfreq[[1]][[i]])
}
min(kmerfreq_rows)

Minus_pre_input = TF_Minus_scalefactors_kmerfreq
for (i in 1:length(Minus_pre_input)) {
  for (p in 1:length(Minus_pre_input[[i]])) {
    Minus_pre_input[[i]][[p]] = Minus_pre_input[[i]][[p]][1:min(kmerfreq_rows)]
  }
}
#Shift each mask by the farthest upstream 3-mer's offset:
itnum = NULL
for (i in 1:40) {
  itnum[i] = 41-i
}

input_correction = data.frame(matrix(ncol = 2))
correction_store = data.frame(matrix(ncol = 2))
for (i in 1:length(itnum)) {
  correction_store[1:itnum[i],1] = rep(1+(i-1), itnum[i]) 
  correction_store[1:itnum[i],2] = rep(137+(i-1), itnum[i]) 
  input_correction = rbind(input_correction, correction_store)
  correction_store = data.frame(matrix(ncol = 2))
}
input_correction = input_correction[-c(1),]
rownames(input_correction) = 1:nrow(input_correction)
#Each mask shifts predictive output by the shift relative to the cut site
#Correct this by assigning the values of each window for shifting
Minus_pred_input = vector(mode = 'list', length = length(Minus_pre_input))
names(Minus_pred_input) <- names(Minus_pre_input)
for (i in 1:length(Minus_pre_input)) {
  for (p in 1:length(Minus_pre_input[[i]])) {
    Minus_pred_input[[i]][[p]] = 
      Minus_pre_input[[i]][[p]][input_correction[p,1]:input_correction[p,2]]
  }
  Minus_pred_input[[i]] = data.frame(Minus_pred_input[[i]])
}
#Assign names
for (i in 1:length(Minus_pred_input)) {
  names(Minus_pred_input[[i]]) = names(Minus_pre_input[[i]])  
}

for (i in 1:length(Minus_pred_input)) {
  names(Minus_pred_input[[i]]) = substr(names(Minus_pred_input[[i]]),
                                       nchar(names(Minus_pred_input[[i]])) - 54, nchar(names(Minus_pred_input[[i]]))-8)
}

###Add x axis to input for future alignment!
for (i in 1:length(Minus_pred_input)) {
  for (p in 1:length(Minus_pred_input[[i]])) {
    Minus_pred_input[[i]][[p]] = as.data.frame(Minus_pred_input[[i]][[p]])
    Minus_pred_input[[i]][[p]][,2] = seq(-73,min(kmerfreq_rows)-73,1)
    colnames(Minus_pred_input[[i]][[p]]) = c('sfkf', 'x')
  }
}

save(Minus_pred_input, file = 'Tn5_TF_Minus_pre_input_spaced6mers.Rdata')
