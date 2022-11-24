source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
###################################################################################
#Benzonase sepcat FASTA files:
uppercasenames <- list('mm39_liver_Benzonase_sepcat.fasta')

#Make sure all FASTA entries are uppercase
uplist <- lapply(uppercasenames, uppercase)  

#Each list object of uplist is the corresponding file:
Benzonase_sepcat = as.data.frame(uplist[1])

rm(uplist, uppercasenames)
#
#Count each nucleotide at each position:
##Benzonase_sepcat
Benzonase_sepcat.transfac = vector("list", length = ceiling(nrow(Benzonase_sepcat)/1000000))
Benzonase_sepcat.split = seqlast(from=1,to=(nrow(Benzonase_sepcat)+1),by=1000000)

for (i in 1:length(Benzonase_sepcat.transfac)) {
  Benzonase_sepcat.transfac[[i]] =
    transfac.func.2(Benzonase_sepcat[Benzonase_sepcat.split[i]:(Benzonase_sepcat.split[i+1]-1),1],31)
}
#Combine all splits of data
transfac.Benzonase_sepcat <- matrix(0,nrow = 31, ncol = 4)
for (i in 1:length(Benzonase_sepcat.transfac)) {
  transfac.Benzonase_sepcat = transfac.Benzonase_sepcat + Benzonase_sepcat.transfac[[i]]
}
#Convert matrix back into DF
transfac.Benzonase_sepcat = cbind(1:nrow(transfac.Benzonase_sepcat), transfac.Benzonase_sepcat)
colnames(transfac.Benzonase_sepcat) = c('P0', 'A', 'C', 'G', 'T')
#Write transfac file for input into weblogo
writeLines(c("ID Benzonase_sepcat_bias",
             "BF Benzonase_sepcat_bias"), 'Benzonase_sepcat_bias.transfac')
write.table(transfac.Benzonase_sepcat, file = "Benzonase_sepcat_bias.transfac",
            append = TRUE, quote=FALSE, row.names =FALSE, col.names = TRUE, sep = '\t')
rm(Benzonase_sepcat, Benzonase_sepcat.transfac)
###################################################################################
###################################################################################
#Cyanase sepcat FASTA file:
uppercasenames <- list('mm39_liver_Cyanase_sepcat.fasta')

#Make sure all FASTA entries are uppercase
uplist <- lapply(uppercasenames, uppercase)  

#Each list object of uplist is the corresponding file:
Cyanase_sepcat = as.data.frame(uplist[1])

rm(uplist, uppercasenames)
##Cyanase_sepcat
Cyanase_sepcat.transfac = vector("list", length = ceiling(nrow(Cyanase_sepcat)/1000000))
Cyanase_sepcat.split = seqlast(from=1,to=(nrow(Cyanase_sepcat)+1),by=1000000)

for (i in 1:length(Cyanase_sepcat.transfac)) {
  Cyanase_sepcat.transfac[[i]] =
    transfac.func.2(Cyanase_sepcat[Cyanase_sepcat.split[i]:(Cyanase_sepcat.split[i+1]-1),1], 31)
}
#Combine all splits of data
transfac.Cyanase_sepcat <- matrix(0,nrow = 31, ncol = 4)
for (i in 1:length(Cyanase_sepcat.transfac)) {
  transfac.Cyanase_sepcat = transfac.Cyanase_sepcat + Cyanase_sepcat.transfac[[i]]
}
#Convert matrix back into DF
transfac.Cyanase_sepcat = cbind(1:nrow(transfac.Cyanase_sepcat), transfac.Cyanase_sepcat)
colnames(transfac.Cyanase_sepcat) = c('P0', 'A', 'C', 'G', 'T')
#Write transfac file for input into weblogo
writeLines(c("ID Cyanase_sepcat_bias",
             "BF Cyanase_sepcat_bias"), 'Cyanase_sepcat_bias.transfac')
write.table(transfac.Cyanase_sepcat, file = "Cyanase_sepcat_bias.transfac",
            append = TRUE, quote=FALSE, row.names =FALSE, col.names = TRUE, sep = '\t')
rm(Cyanase_sepcat.transfac, Cyanase_sepcat)
###################################################################################
###################################################################################
#MNase sepcat FASTA file:
uppercasenames <- list('mm39_liver_MNase_sepcat.fasta')

#Make sure all FASTA entries are uppercase
uplist <- lapply(uppercasenames, uppercase)  

#Each list object of uplist is the corresponding file:
MNase_sepcat = as.data.frame(uplist[1])

rm(uplist, uppercasenames)
##MNase_sepcat
MNase_sepcat.transfac = vector("list", length = ceiling(nrow(MNase_sepcat)/1000000))
MNase_sepcat.split = seqlast(from=1,to=(nrow(MNase_sepcat)+1),by=1000000)

for (i in 1:length(MNase_sepcat.transfac)) {
  MNase_sepcat.transfac[[i]] =
    transfac.func.2(MNase_sepcat[MNase_sepcat.split[i]:(MNase_sepcat.split[i+1]-1),1],
                    'MNase_sepcat_bias', 31, FALSE)
}
#Combine all splits of data
transfac.MNase_sepcat <- matrix(0,nrow = 31, ncol = 4)
for (i in 1:length(MNase_sepcat.transfac)) {
  transfac.MNase_sepcat = transfac.MNase_sepcat + MNase_sepcat.transfac[[i]]
}
#Convert matrix back into DF
transfac.MNase_sepcat = cbind(1:nrow(transfac.MNase_sepcat), transfac.MNase_sepcat)
colnames(transfac.MNase_sepcat) = c('P0', 'A', 'C', 'G', 'T')
#Write transfac file for input into weblogo
writeLines(c("ID MNase_sepcat_bias",
             "BF MNase_sepcat_bias"), 'MNase_sepcat_bias.transfac')
write.table(transfac.MNase_sepcat, file = "MNase_sepcat_bias.transfac",
            append = TRUE, quote=FALSE, row.names =FALSE, col.names = TRUE, sep = '\t')
rm(MNase_sepcat.transfac, MNase_sepcat)
###################################################################################
###################################################################################
#DNase sepcat FASTA file:
uppercasenames <- list('DNase_Naked_unscaled_sepcat.fasta')

#Make sure all FASTA entries are uppercase
uplist <- lapply(uppercasenames, uppercase)  

DNase_sepcat = as.data.frame(uplist[1])

rm(uplist, uppercasenames)
##DNase_sepcat
DNase_sepcat.transfac = vector("list", length = ceiling(nrow(DNase_sepcat)/1000000))
DNase_sepcat.split = seqlast(from=1,to=(nrow(DNase_sepcat)+1),by=1000000)

for (i in 1:length(DNase_sepcat.transfac)) {
  DNase_sepcat.transfac[[i]] =
    transfac.func.2(DNase_sepcat[DNase_sepcat.split[i]:(DNase_sepcat.split[i+1]-1),1],
                    'DNase_sepcat_bias', 31, FALSE)
}
#Combine all splits of data
transfac.DNase_sepcat <- matrix(0,nrow = 31, ncol = 4)
for (i in 1:length(DNase_sepcat.transfac)) {
  transfac.DNase_sepcat = transfac.DNase_sepcat + DNase_sepcat.transfac[[i]]
}
#Convert matrix back into DF
transfac.DNase_sepcat = cbind(1:nrow(transfac.DNase_sepcat), transfac.DNase_sepcat)
colnames(transfac.DNase_sepcat) = c('P0', 'A', 'C', 'G', 'T')
#Write transfac file for input into weblogo
writeLines(c("ID DNase_sepcat_bias",
             "BF DNase_sepcat_bias"), 'DNase_sepcat_bias.transfac')
write.table(transfac.DNase_sepcat, file = "DNase_sepcat_bias.transfac",
            append = TRUE, quote=FALSE, row.names =FALSE, col.names = TRUE, sep = '\t')
rm(DNase_sepcat.transfac,DNase_sepcat)
###################################################################################
###################################################################################
#Tn5 sepcat FASTA file:
uppercasenames <- list('C1_gDNA_rep1_sepcat.fasta')
#Make sure all FASTA entries are uppercase
uplist <- lapply(uppercasenames, uppercase)  
Tn5_sepcat = as.data.frame(uplist[1])

rm(uplist, uppercasenames)
##Tn5_sepcat
Tn5_sepcat.transfac = vector("list", length = ceiling(nrow(Tn5_sepcat)/1000000))
Tn5_sepcat.split = seqlast(from=1,to=(nrow(Tn5_sepcat)+1),by=1000000)

for (i in 1:length(Tn5_sepcat.transfac)) {
  Tn5_sepcat.transfac[[i]] =
    transfac.func.2(Tn5_sepcat[Tn5_sepcat.split[i]:(Tn5_sepcat.split[i+1]-1),1],
                    'Tn5_sepcat_bias', 31, FALSE)
}
#Combine all splits of data
transfac.Tn5_sepcat <- matrix(0,nrow = 31, ncol = 4)
for (i in 1:length(Tn5_sepcat.transfac)) {
  transfac.Tn5_sepcat = transfac.Tn5_sepcat + Tn5_sepcat.transfac[[i]]
}
#Convert matrix back into DF
transfac.Tn5_sepcat = cbind(1:nrow(transfac.Tn5_sepcat), transfac.Tn5_sepcat)
colnames(transfac.Tn5_sepcat) = c('P0', 'A', 'C', 'G', 'T')
#Write transfac file for input into weblogo
writeLines(c("ID Tn5_sepcat_bias",
             "BF Tn5_sepcat_bias"), 'Tn5_sepcat_bias.transfac')
write.table(transfac.Tn5_sepcat, file = "Tn5_sepcat_bias.transfac",
            append = TRUE, quote=FALSE, row.names =FALSE, col.names = TRUE, sep = '\t')
###################################################################################