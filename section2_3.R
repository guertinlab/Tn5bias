source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
library(ggseqlogo)

pswm.func.2 <- function(x.ligation, out = 'outfilename', posnum, rc = FALSE) {
   col.matrix = matrix()
   for (g in 1:posnum){
    itnum = lapply(strsplit(as.character(x.ligation), ''), "[", g)
    if (g == 1) {
    col.matrix = itnum
    } else {
    col.matrix = cbind(col.matrix, itnum)
    }
    }  
  
  a.nuc = sapply(1:posnum, function(x) sum(col.matrix[,x] == "A"))
  t.nuc = sapply(1:posnum, function(x) sum(col.matrix[,x] == "T"))
  c.nuc = sapply(1:posnum, function(x) sum(col.matrix[,x] == "C"))
  g.nuc = sapply(1:posnum, function(x) sum(col.matrix[,x] == "G"))
  
  pswm = cbind(a.nuc, c.nuc, g.nuc, t.nuc)
  print(pswm)
  outfile = file(paste0(out, '.txt'))
  on.exit(close(outfile))
  writeLines(c("MEME version 4", "ALPHABET= ACGT", "strands: + -", " ", 
               "Background letter frequencies (from uniform background):", 
               "A 0.30000 C 0.20000 G 0.20000 T 0.30000", paste("MOTIF", out), " ",
               "letter-probability matrix: alength= 4 w= 20"), outfile)
  pswm = pswm/rowSums(pswm)
  if (rc == "TRUE"){
    pswm<- pswm[nrow(pswm):1,ncol(pswm):1]
  } else {}
      
  
  write.table(pswm, file = paste0(out, '.txt'), append = TRUE, quote=FALSE, row.names =FALSE, col.names = FALSE)
  return(pswm)
}


plot.seqlogo.func <- function(x, outfile = "ATAC-kmer_optimization_all_test.pdf") {
  w =  0.663 + (ncol(x) + 1)*0.018 + (ncol(x)+2)* .336
  pdf(outfile, useDingbats=FALSE, width=w, height=2.695)
  print(ggseqlogo(x,  facet = "wrap", font = 'helvetica_bold'))
  dev.off()
}

dir = paste0('Your Working Directory Here')
setwd(dir)

#
#Input filenames (from wd) into uppercasenames in this order: pe1.minus, pe1.plus, pe2.minus, pe2.plus
uppercasenames <- list('pe1.minus', 'pe1.plus', 'pe2.minus', 'pe2.plus')

uppercase <- function(tableinput){ 
                      upperdf <- read.table(tableinput, comment.char = '>')
                      upperdf[,1] = as.character(upperdf[,1])
                      upperdf = data.frame(lapply(upperdf, function(v) {
                        if (is.character(v)) return(toupper(v))
                        else return(v)}))
                        }

uplist <- lapply(uppercasenames, uppercase)  
  
pe1.minusATAC = as.data.frame(uplist[1])
pe1.plusATAC = as.data.frame(uplist[2])
pe2.minusATAC = as.data.frame(uplist[3])
pe2.plusATAC = as.data.frame(uplist[4])

pswm.pe1.minus = pswm.func.2(pe1.minusATAC[,1], 'ATAC_bias_pe1_minus', 21, FALSE)
pswm.pe1.plus = pswm.func.2(pe1.plusATAC[,1], 'ATAC_bias_pe1_plus', 21, FALSE)
pswm.pe2.minus = pswm.func.2(pe2.minusATAC[,1], 'ATAC_bias_pe2_minus', 21, FALSE)
pswm.pe2.plus = pswm.func.2(pe2.plusATAC[,1], 'ATAC_bias_pe2_plus', 21, FALSE)

save(pswm.pe1.minus, pswm.pe1.plus, pswm.pe2.minus, pswm.pe2.plus, file ='pswm.Rdata')


pswm.pe1.minus.trans = t(pswm.pe1.minus)
rownames(pswm.pe1.minus.trans) = c('A', 'C', 'G', 'T')

pswm.pe1.plus.trans = t(pswm.pe1.plus)
rownames(pswm.pe1.plus.trans) = c('A', 'C', 'G', 'T')

pswm.pe2.minus.trans = t(pswm.pe2.minus)
rownames(pswm.pe2.minus.trans) = c('A', 'C', 'G', 'T')

pswm.pe2.plus.trans = t(pswm.pe2.plus)
rownames(pswm.pe2.plus.trans) = c('A', 'C', 'G', 'T')

plot.seqlogo.func(pswm.pe1.minus.trans, outfile='pswm_pe1_minus_trans.pdf')
plot.seqlogo.func(pswm.pe1.plus.trans, outfile='pswm_pe1_plus_trans.pdf')

plot.seqlogo.func(pswm.pe2.minus.trans, outfile='pswm_pe2_minus_trans.pdf')
plot.seqlogo.func(pswm.pe2.plus.trans, outfile='pswm_pe2_plus_trans.pdf')
