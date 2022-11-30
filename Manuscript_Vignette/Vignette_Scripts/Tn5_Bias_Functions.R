#Calculate scale factors from a seqOutBias table file
scalefactor.func <- function(factortable){
    require(data.table)
    datatable = fread(factortable, skip = 1)
    datatable[, pluscountpercent := datatable[,3]/sum(datatable[,3])]
    datatable[, plusobspercent := datatable[,5]/sum(datatable[,5])]
    datatable[, plusscalefact := round(datatable[,pluscountpercent]/datatable[,plusobspercent], 6)]
    datatable[, minuscountpercent := datatable[,4]/sum(datatable[,4])]
    datatable[, minusobspercent := datatable[,6]/sum(datatable[,6])]
    datatable[, minusscalefact := round(datatable[,minuscountpercent]/datatable[,minusobspercent], 6)]
}
#

#FIMO.to.Bed takes input in FIMO format (1-based start and end) and converts 
#it to BED format (0-based start, 1-based end), also increases window by amount specified.
FIMO.to.BED <- function(fimofile, window = 0) {
    require(data.table)
    #load in FIMO file
    FIMOdf <- fread(fimofile, sep = '\t', header = TRUE)
    #shift start position by -1 to make it 0-based
    FIMOdf[,4] = FIMOdf[,4] - 1
    #make location column
    FIMOdf[,location := paste(sequence_name, ':', start, '-', stop, sep = '')]
    #Rearrange columns for BED format
    beddf <- FIMOdf[, c(3,4,5,2,11,6)]
    colnames(beddf) <- c('chr', 'start', 'end', 'gene', 'location', 'strand')
    #CB (central base) determines the central base which the window is based on
    #also used as the 0 position on x-axis
    #even motifs will == 0, odd will == 1
    if ((beddf[1,3] - beddf[1,2])%%2 == 0) {
        CB = data.table(beddf)
        #convert start positions to 1-base for calculating CB
        CB[,start := start + 1]
        #Calculate median ceiling for positive motifs, and median floor value
        # for minus motifs then use this value to determine
        #start and end using upstream/downstream
        CB[,median := apply(CB[,2:3],1, median, na.rm = TRUE)][]
        setkey(CB,strand)
        CB[c("-"),median := floor(median)] 
        CB[c("+"),median := ceiling(median)]
        CB[,start := median - window]
        CB[,end := median + window]
        #Convert start position back to 0-base
        CB[,start := start - 1]
        CB[,median:=NULL]
    } else {
        CB = data.table(beddf)
        #convert start positions to 1-base for calculating median
        CB[,start := start + 1]
        #Calculate median for each row, then use this value to determine
        #start and end using upstream/downstream
        CB[,median := apply(CB[,2:3],1, median, na.rm = TRUE)][]
        CB[,start := median - window]
        CB[,end := median + window]
        #Convert start position back to 0-base
        CB[,start := start - 1]
        CB[,median:=NULL]
    }

    return(CB)
}

#Split up a list of sequences into mers = mermask for each mappable position in the sequence
mer.positions <- function(input, mermask = 'NNNXXXXXNNN') {
    seqdf_mer <- vector(mode = "list", length = nchar(input[1]))
    if (grepl('X', mermask) == TRUE) {
    mermask_spaces = unlist(strsplit(mermask, 'X'))
    mer1 = nchar(mermask_spaces)[1]
    unmasked = length(which(mermask_spaces == '')) + 1
    mer2 = nchar(mermask_spaces)[length(mermask_spaces)]
    for (i in 1:nchar(input[1])) {
        seqdf_store = data.table()
        seqdf_store[, alltogether := paste(substr(input[1:length(input)], i, i + mer1 - 1),
                                           paste(rep('X', unmasked), collapse = ""),
                                           substr(input[1:length(input)],
                                           i + mer1 + unmasked, i + mer1 + unmasked + mer2 - 1), sep = '')]
        seqdf_mer[[i]] <- seqdf_store$alltogether 
    }
    mersize = nchar(mermask)
    names(seqdf_mer) <- paste0("Position_", seq_along(1:nchar(input[1])))
    seqdf_mer <- seqdf_mer[1:(length(seqdf_mer)-mersize+1)]}
    else {
      seqdf_store <- NULL
      mersize = nchar(mermask)
      for (i in 1:nchar(input[1])) {
        seqdf_store <- substr(input[1:length(input)], i, i+ (mersize-1))
        seqdf_mer[[i]] <- seqdf_store 
        names(seqdf_mer) <- paste0("Position_", seq_along(1:nchar(input[1])))
      }
      seqdf_mer <- seqdf_mer[1:(length(seqdf_mer)-mersize+1)]
    }
    return(seqdf_mer)
}


#Get mercounts for each possible mer from the output of mer.positions - only for 1 sequence
mer.counts = function(sequence_mer_positions, mermask = 'NNNN') {
  #Create all possible combinations of kmer
  mermask = unlist(strsplit(mermask, split = ''))
  mertable = expand.grid(rep(list(c('A','C','G','T')), length(which(mermask == 'N'))))
  if (length(which(mermask == 'X')) > 0) {
    Xrows = matrix(nrow = nrow(mertable), ncol = length(which(mermask == 'X')))
    Xrows[,c(1:length(which(mermask == 'X')))] = rep('X', nrow(mertable))
    mertable = cbind(mertable[,1:(which(mermask=='X')[1]-1)],
                     Xrows, mertable[,(which(mermask=='X')[1]:ncol(mertable))])
  } else {}
  mertable = data.frame(apply(mertable, 1 , paste, collapse = ""))
  mercount = as.character(sequence_mer_positions)
  merstore = NULL
  #
  for (i in 1:nrow(mertable)) {
    merstore = length(which(mercount %in% mertable[i,1]))
    if (length(merstore) == 0) {
      merstore = 0
    } else {}
    mertable[i,2] = merstore
  }
  colnames(mertable) = c('kmer', 'count')
  return(mertable)
}


#Get position frequencies for each possible mer from the output of mer.positions
position.frequencies = function(sequence_mer_positions, mermask = 'NNNNNNN') {
  #Create all possible combinations of kmer
  mermask = unlist(strsplit(mermask, split = ''))
  mertable = expand.grid(rep(list(c('A','C','G','T')), length(which(mermask == 'N'))))
  if (length(which(mermask == 'X')) > 0) {
  Xrows = matrix(nrow = nrow(mertable), ncol = length(which(mermask == 'X')))
  Xrows[,c(1:length(which(mermask == 'X')))] = rep('X', nrow(mertable))
  mertable = cbind(mertable[,1:(which(mermask=='X')[1]-1)],
                   Xrows, mertable[,(which(mermask=='X')[1]:ncol(mertable))])
  } else {}
  mertable = data.table(apply(mertable, 1 , paste, collapse = ""))
  posfreq_poslist = vector(mode = "list", length = length(sequence_mer_positions))
  #
  for (p in 1:length(sequence_mer_positions)) {
    posfreq_posstore = data.table(mertable$V1)
    postable = table(sequence_mer_positions[[p]])
    if (length(which(!names(postable) %in% mertable$V1)) > 0) {
    postable = postable[-c(which(!names(postable) %in% mertable$V1))]
    } else {}
    missing_mers = as.character(mertable$V1[which(!mertable$V1 %in% names(postable))])
    append_mers = rep(0, length(missing_mers))
    names(append_mers) = missing_mers
    postable = c(postable, append_mers)
    postable = postable/length(sequence_mer_positions[[1]])
    posfreq_poslist[[p]] <- postable
    
  }
  names(posfreq_poslist) <- paste0("Position_", seq_along(1:length(posfreq_poslist)))
  return(posfreq_poslist)
}
##Multiply 1/scale_factor by its corresponding kmer frequency for each position in
##the composite
scalefactor.by.kmerfrequency <- function(scalefactors, kmerfrequency, tfname) {
    all_scalefactors = vector(mode = "list", length = length(scalefactors))
    names(all_scalefactors) <- names(scalefactors)
   for (j in 1:length(scalefactors)) {
       scalefactor = scalefactors[[j]]
       scalefactor = scalefactor[order(V2)]
       scaledkfreq = vector(mode = "list", length = length(kmerfrequency))
       names(scaledkfreq) = names(kmerfrequency)
    for (i in 1:length(kmerfrequency)) {
        scaledkfreq[[i]] = data.table(stack(kmerfrequency[[i]]))
        colnames(scaledkfreq[[i]]) = c('kmerfreq', 'kmer')
        scaledkfreq[[i]]$kmer = as.character(scaledkfreq[[i]]$kmer)
        scaledkfreq[[i]] = scaledkfreq[[i]][order(kmer)]
        scaledkfreq[[i]][, mask := scalefactor[,2]]  
        scaledkfreq[[i]][, plusscalefreq := (1/(scalefactor[,9])) * kmerfreq ]
        scaledkfreq[[i]][, minusscalefreq := (1/(scalefactor[,12])) * kmerfreq ]
    
        }
       saveRDS(scaledkfreq, file = paste(tfname, '_',
              gsub('_scale_factors.txt', '', names(scalefactors[j])), '.rds', sep = ''))
       #all_scalefactors[[j]] = scaledkfreq 
       }
      #return(all_scalefactors)
}

##Take sum of each position's scaledkmerfreq (1/scale_factor * kmer frequency) and use to create 
##data.table of inputs for each mask for each TF motif
scaledkmerfreq.sum.plus <- function(input) {
    output = data.table(matrix(nrow = length(input)))
    for (i in 1:length(input)) {

            output[i, scaledkmerfreqsumplus := sum(input[[i]][,4])]

    }
    output[,1] <- NULL
    rownames(output) = names(input)
    return(output)
}

scaledkmerfreq.sum.minus <- function(input) {
  output = data.table(matrix(nrow = length(input)))
  for (i in 1:length(input)) {
    
    output[i, scaledkmerfreqsumminus := sum(input[[i]][,5])]
    
  }
  output[,1] <- NULL
  rownames(output) = names(input)
  return(output)
}



################################################################################
################################################################################
#Combine values from mask positions to plot redudant importance values for each position
redundantpos = function(input, labelnum = 2, impnum = 3) {
  input[,labelnum] = as.character(input[,labelnum])
  input[,labelnum] = gsub('C', '', input[,labelnum])
  posimp = vector(mode = 'list', length = nchar(input[1,labelnum]))
  for (i in 1:nchar(input[1,labelnum])) {
    for (p in 1:nrow(input)) {
      posimp[[i]][p] = fifelse(substr(input[p,labelnum],i,i) == 'N', input[p,impnum],0)
    }
  }
  output = data.table(matrix(nrow = length(posimp), ncol = 1))
  for (i in 1:length(posimp)) {
    output[i, possum := sum(posimp[[i]])]
  }
  output = output[,-c(1)]
  output[, position := paste(rep('X',36,sep = ''), collapse = '')]
  
  for (i in 1:nrow(output)) {
    output[i,]$position = paste(substr(output[i, position],1,i-1),'N', 
                                substr(output[i, position],i+1,nchar(output[i, position])), sep = '')
    output[i,]$position = paste(substr(output[i, position],1,(nchar(output[i, position])/2)), 'C', 
                                paste(substr(output[i, position],
                                (nchar(output[i, position])/2+1),nchar(output[i, position]))), sep = '')
  }
  return(output)
}


#
#BED.query.bigWig uses the bed coordinates produced by FIMO.to.BED to
#query a supplied bigWig file for the read counts at the specified positions
#In order to create a window around the base of interest, a central base (CB)
#is determined based on whether or not the motif of interest is odd or even
#additionally, if ATAC=FALSE, minus-aligned motifs are shifted +1 
#relative to the plus strand to account for DNase etc. cutting between bases
#finally, the mean of each position's read counts is calculated for plotting
BED.query.bigWig <- function(beddf, bwPlus, bwMinus, upstream = 10,
                             downstream = 10,
                             factor = '', group = '', ATAC = TRUE) {
  require(bigWig)
  require(data.table)
  require(matrixStats)
  step = 1
  #load bigWigs
  bw.plus = load.bigWig(bwPlus)
  bw.minus = load.bigWig(bwMinus)
  #CB (central base) determines the central base which the window is based on
  #also used as the 0 position on x-axis
  #even motifs will == 0, odd will == 1
  if ((beddf[1,3] - beddf[1,2])%%2 == 0) {
    CB = data.table(beddf)
    #convert start positions to 1-base for calculating CB
    CB[,start := start + 1]
    #Calculate median ceiling for positive motifs, and median floor value
    # for minus motifs then use this value to determine
    #start and end using upstream/downstream
    CB[,median := apply(CB[,2:3],1, median, na.rm = TRUE)][]
    setkey(CB,strand)
    CB[c("-"),median := floor(median)] 
    CB[c("+"),median := ceiling(median)]
    CB[,start := median - upstream]
    CB[,end := median + downstream]
    #Convert start position back to 0-base
    CB[,start := start - 1]
    CB[,median:=NULL]
  } else {
    CB = data.table(beddf)
    #convert start positions to 1-base for calculating median
    CB[,start := start + 1]
    #Calculate median for each row, then use this value to determine
    #start and end using upstream/downstream
    CB[,median := apply(CB[,2:3],1, median, na.rm = TRUE)][]
    CB[,start := median - upstream]
    CB[,end := median + downstream]
    #Convert start position back to 0-base
    CB[,start := start - 1]
    CB[,median:=NULL]
  }
  #Shift minus-aligned motifs by +1 for DNase etc.
  if (ATAC == FALSE) {
    setkey(CB,strand)
    CB[c("-"),start := start + 1] 
    CB[c("-"),end := end + 1] 
  } else {}
  #extract read count information from bigWigs, using CB-derived
  #window as coordinates
  tss.matrix = bed6.step.bpQuery.bigWig(bw.plus, bw.minus , CB,
                                        step = 1, as.matrix=TRUE, follow.strand=TRUE)
  tss_correction = sum(tss.matrix)/(nrow(tss.matrix)*ncol(tss.matrix))
  #Determine beginning and end of plot
  #if not ATAC, peaks are between bases, so offset by -1 and 
  #add 0.5 to align 'cuts' on 0 at -0.5
  if (ATAC == FALSE) {
    coordin.start = (-upstream - 0.5 )
    coordin.end = (downstream - 0.5 )
  } else {
    #For tn5, add negative upstream value and downstream value
    coordin.start = (-upstream)
    coordin.end = (downstream)
  }
  #create coordinates for aligning read counts on composite profile
  #also take means of each column from tss.matrix and label each row with
  #factor column
  composite.lattice = data.table(seq(coordin.start, coordin.end, by = step),
                                 colMeans(tss.matrix),
                                 factor, group, tss_correction, 
                                 stringsAsFactors = FALSE)
  colnames(composite.lattice) = c('x', 'est', 'factor', 'group', 'correction')
  composite.lattice$x = as.numeric(composite.lattice$x)
  unload.bigWig(bw.plus)
  unload.bigWig(bw.minus)
  return(composite.lattice)
}
#

##Plot.composites takes the composite.lattice object and creates a plot for 
#this data, while also allowing a mapping of the original motif width onto the
#plot by specifying a TFlen value (the default is 10). 
plot.composites <- function(dat, x_axis_range = min(dat$x):max(dat$x), ylabel = '',
                            pdf_name = 'PLEASE_SET_FILE_NAME', y_axis_range = min(dat$est):max(dat$est),
                            xlabel = '', striplabel = TRUE, legend = TRUE,
                            motifline = FALSE, Motiflen = 10, nYticks = 3,
                            figwidth = 2.5, figheight=3, hline = FALSE, nYaxisdigits = 3,
                            indexlist = NULL, layoutgrid = NULL, y_axis = FALSE, linethick = 2,
                            X_axis_ticks = seq(-20,20,10), hline_val = 0, Y_axis_ticks = seq(0,0.02,0.01),
                    col.lines = c("#0000FF", "#FF0000", "#00000090", 
                                  rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                  rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), 
                                  rgb(1/2,1/2,0,1/2)), 
                                  fill.poly = c(rgb(0,0,1,1/4), 
                                  rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4), 
                                  rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    require(lattice)
    dat = dat[which(dat$x <= max(x_axis_range)+0.5 & dat$x >= min(x_axis_range) -0.5),]
    
    pdf(paste(pdf_name, '.pdf', sep = ''), width= figwidth, height= figheight) 
    print(xyplot(est ~ x|factor, group = group, data = dat, strip = striplabel,
                 type = 'l', as.table = TRUE,
                 scales=list(x=list(at= X_axis_ticks,cex=1.1, relation = "free", axs ="i", rot = 45), 
                             y =list(at=as.numeric(format(Y_axis_ticks, digits = nYaxisdigits)),
                               cex=1.1, relation = "free", rot = 0)),
                 col = col.lines,
                 auto.key = if (legend == TRUE) 
                     {list(points=F, lines=T, cex=0.8)} else{},
                 par.settings = list(strip.background=list(col="#00000000"),
                                     strip.border = list(col = 'transparent'),
                                     superpose.symbol = list(pch = c(16),
                                                             col=col.lines, 
                                                             cex =0.5), 
                                     superpose.line = list(col = col.lines, 
                                                           lwd=c(2), 
                                                           lty = c(1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=1.2, font=1, col='black'),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 lwd=linethick,
                 ylim = if (y_axis == TRUE){c(min(y_axis_range), max(y_axis_range))} else{},
                 xlim = c(min(x_axis_range), max(x_axis_range)),
                 ylab = list(label = paste(ylabel), cex =1.3),
                 xlab = list(label = paste(xlabel), cex =1.3),
                 index.cond = indexlist,
                 layout = layoutgrid,
                 panel = function(x, y, ...) {
                   panel.xyplot(x, y, ...)
                   if (hline == TRUE)
                   {panel.abline(h = hline_val, lty =2, lwd = 2, col = '#FF0000')}else{}
                   level = dimnames(trellis.last.object())[["factor"]][packet.number()]
                   if (motifline == TRUE) 
                   {panel.abline(v = Motiflen[rownames(Motiflen)==level,]/2, lty = 2, col = "red")} else{}
                   if (motifline == TRUE) 
                   {panel.abline(v = -Motiflen[rownames(Motiflen)==level,]/2, lty = 2, col = "red")} else{}
                 }
    ))
    dev.off()
}
#
##Plots seqlogos
plot.seqlogo.func <- function(x, outfile = "PLEASE_SET_FILE_NAME.pdf") {
  require(ggseqlogo)
  w =  0.663 + (ncol(x) + 1)*0.018 + (ncol(x)+2)* .336
  pdf(outfile, useDingbats=FALSE, width=w, height=2.695)
  print(ggseqlogo(x,  facet = "wrap", font = 'helvetica_bold'))
  dev.off()
}


#Count all instances of each nucleotide at each position, data frame to be counted is
#x.ligation, and posnum is the size of the window counted
transfac.func.2 <- function(x.ligation, posnum) {
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
  
  transfac = cbind(a.nuc, c.nuc, g.nuc, t.nuc)
  print(transfac)
  return(transfac)
}

#Split up data sets into groups of desired size to reduce RAM load
seqlast <- function (from, to, by) 
{
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

#Make all instaces in a table uppercase
uppercase <- function(tableinput){ 
  upperdf <- read.table(tableinput, comment.char = '>')
  upperdf[,1] = as.character(upperdf[,1])
  upperdf = data.frame(lapply(upperdf, function(v) {
    if (is.character(v)) return(toupper(v))
    else return(v)}))
}

#Plots seqlogos of JASPAR PSWMs
plot.seqlogo.meme <-function(infile, outfile = "memefunc_test.pdf") {
  require(ggseqlogo)  
  wd <- getwd()
  memefile <- paste0(wd,"/",infile)
  minmeme <- read.csv(memefile, skip = 11, header = FALSE, sep = '')
  minmeme <- minmeme[1:nrow(minmeme)-1,]
  meme.trans <- t(minmeme)
  class(meme.trans) <- "numeric"
  rownames(meme.trans) = c('A', 'C', 'G', 'T')
  w =  0.663 + (ncol(meme.trans) + 1)*0.018 + (ncol(meme.trans)+2)* .336
  pdf(outfile, useDingbats=FALSE, width=w, height=2.695)
  print(ggseqlogo(meme.trans,  facet = "wrap", font = 'helvetica_bold'))
  dev.off()
}

#Import rule ensemble formula in usable form from output .txt file                 
import.rules = function(input) {
  rulesdf = data.table(read.table(input, sep = ',', header = TRUE))
  lincoeff = NULL
  rulescoeff = NULL
  for (i in 1:nrow(rulesdf)) {
    lincoeff[i] = fifelse(!grepl(">|<|=",rulesdf[i, description]), 
                          paste('(',rulesdf[i, description], '*',
                                rulesdf[i, coefficient], ')', sep = ''),paste(''))}
  lincoeff = lincoeff[!(lincoeff=='')]
  linpaste = capture.output(cat(lincoeff, sep = ' + '))
  for (i in 1:nrow(rulesdf)) {
    rulescoeff[i] = fifelse(grepl(">|<|=",rulesdf[i, description]),
                            paste('fifelse(', rulesdf[i, description], ', ', 
                                  rulesdf[i, coefficient], ', 0)', sep = ''), '')}
  rulescoeff = rulescoeff[!(rulescoeff=='')]
  rulespaste = capture.output(cat(rulescoeff, sep = ' + '))
  output = paste(linpaste, rulespaste, sep = ' + ')
  return(output)
}                 

#Add violin plot on top of other lattice objects
panel.violin.hack <-
  function (x, y, box.ratio = 1, box.width = box.ratio/(1 + box.ratio),
            horizontal = TRUE, alpha = plot.polygon$alpha, border =  
              plot.polygon$border,
            lty = plot.polygon$lty, lwd = plot.polygon$lwd, col = plot.polygon 
            $col,
            varwidth = FALSE, bw = NULL, adjust = NULL, kernel = NULL,
            window = NULL, width = NULL, n = 50, from = NULL, to = NULL,
            cut = NULL, na.rm = TRUE, ...){
    if (all(is.na(x) | is.na(y)))
      return()
    x <- as.numeric(x)
    y <- as.numeric(y)
    plot.polygon <- trellis.par.get("plot.polygon")
    darg <- list()
    darg$bw <- bw
    darg$adjust <- adjust
    darg$kernel <- kernel
    darg$window <- window
    darg$width <- width
    darg$n <- n
    darg$from <- from
    darg$to <- to
    darg$cut <- cut
    darg$na.rm <- na.rm
    my.density <- function(x) {
      ans <- try(do.call("density", c(list(x = x), darg)),
                 silent = TRUE)
      if (inherits(ans, "try-error"))
        list(x = rep(x[1], 3), y = c(0, 1, 0))
      else ans
    }
    numeric.list <- if (horizontal)
      split(x, factor(y))
    else split(y, factor(x))
    levels.fos <- as.numeric(names(numeric.list))
    d.list <- lapply(numeric.list, my.density)
    dx.list <- lapply(d.list, "[[", "x")
    dy.list <- lapply(d.list, "[[", "y")
    max.d <- sapply(dy.list, max)
    if (varwidth)
      max.d[] <- max(max.d)
    xscale <- current.panel.limits()$xlim
    yscale <- current.panel.limits()$ylim
    height <- box.width
    if (horizontal) {
      for (i in seq_along(levels.fos)) {
        if (is.finite(max.d[i])) {
          pushViewport(viewport(y = unit(levels.fos[i],
                                         "native"), height = unit(height, "native"),
                                yscale = c(max.d[i] * c(-1, 1)), xscale = xscale))
          grid.polygon(x = c(dx.list[[i]], rev(dx.list[[i]])),
                       y = c(dy.list[[i]], -rev(dy.list[[i]])),  
                       default.units = "native",
                       # this is the point at which the index is added
                       gp = gpar(fill = col[i], col = border, lty = lty,
                                 lwd = lwd, alpha = alpha))
          popViewport()
        }
      }
    }
    else {
      for (i in seq_along(levels.fos)) {
        if (is.finite(max.d[i])) {
          pushViewport(viewport(x = unit(levels.fos[i],
                                         "native"), width = unit(height, "native"),
                                xscale = c(max.d[i] * c(-1, 1)), yscale = yscale))
          grid.polygon(y = c(dx.list[[i]], rev(dx.list[[i]])),
                       x = c(dy.list[[i]], -rev(dy.list[[i]])),  
                       default.units = "native",
                       # this is the point at which the index is added
                       gp = gpar(fill = col[i], col = border, lty = lty,
                                 lwd = lwd, alpha = alpha))
          popViewport()
        }
      }
    }
    invisible()
  }
                 
            
