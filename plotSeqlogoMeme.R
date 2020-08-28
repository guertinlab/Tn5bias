library(ggseqlogo)

plot.seqlogo.meme <-function(infile, outfile = "memefunc_test.pdf")
{
  wd <- getwd()
  memefile <- paste0(wd,"/",infile)
  minmeme <- read.csv(memefile, skip = 12, header = FALSE, sep = '')
  minmeme <- minmeme[1:nrow(minmeme)-1,]
  meme.trans <- t(minmeme)
  class(meme.trans) <- "numeric"
  rownames(meme.trans) = c('A', 'C', 'G', 'T')
  w =  0.663 + (ncol(meme.trans) + 1)*0.018 + (ncol(meme.trans)+2)* .336
  pdf(outfile, useDingbats=FALSE, width=w, height=2.695)
  print(ggseqlogo(meme.trans,  facet = "wrap", font = 'helvetica_bold'))
  dev.off()
}


plot.seqlogo.meme('XXXXXX_minimal_meme.txt', outfile = "XXXXXX_seqlogo.pdf")
