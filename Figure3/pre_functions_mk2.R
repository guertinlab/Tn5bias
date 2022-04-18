
#
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

#Split up a list of sequences into mers = mersize for each mappable position in the sequence
mer.positions <- function(input, mersize = 5) {
    seqdf_mer <- vector(mode = "list", length = nchar(input[1]))
    seqdf_store <- NULL  
    for (i in 1:nchar(input[1])) {
        seqdf_store <- substr(input[1:length(input)], i, i+ (mersize-1))
        seqdf_mer[[i]] <- seqdf_store 
        names(seqdf_mer) <- paste0("Position_", seq_along(1:nchar(input[1])))
    }
    seqdf_mer <- seqdf_mer[1:(length(seqdf_mer)-mersize+1)]
    return(seqdf_mer)
}

#Get position frequencies for each possible mer from the output of mer.positions
position.frequencies <- function(sequence_mer_positions) {
    #Create all possible combinations of kmer
    mertable = expand.grid(rep(list(c('A','C','G','T')), nchar(sequence_mer_positions[[1]][1])))
    mertable = data.table(apply(mertable, 1 , paste, collapse = ""))
    posfreq_posstore = NULL
    posfreq_poslist = vector(mode = "list", length = length(sequence_mer_positions))
    #
    for (i in 1:length(sequence_mer_positions)) {
        for (p in 1:length(mertable$V1)) {
            posfreq_posstore[p] = length(which(sequence_mer_positions[[i]] == mertable[[1]][p]))/length(sequence_mer_positions[[i]])
        }
        posfreq_poslist[[i]] <- posfreq_posstore
        names(posfreq_poslist[[i]]) <- mertable$V1
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
       scaledkfreq = vector(mode = "list", length = length(kmerfrequency))
       names(scaledkfreq) = names(kmerfrequency)
    for (i in 1:length(kmerfrequency)) {
        scaledkfreq[[i]] = data.table(stack(kmerfrequency[[i]]))
        colnames(scaledkfreq[[i]]) <- c('kmerfreq', 'kmer')
        scaledkfreq[[i]]$kmer <- as.character(scaledkfreq[[i]]$kmer)
        for (p in 1:length(kmerfrequency[[i]])) {
            scaledkfreq[[i]][p, mask := scalefactor[grep(scaledkfreq[[i]]$kmer[p], scalefactor$V2),2]]  
            scaledkfreq[[i]][p, plusscalefreq := (1/(scalefactor[grep(scaledkfreq[[i]]$kmer[p], scalefactor$V2),9])) * kmerfreq ]
            scaledkfreq[[i]][p, minusscalefreq := (1/(scalefactor[grep(scaledkfreq[[i]]$kmer[p], scalefactor$V2),12])) * kmerfreq ]
            }
        }
       saveRDS(scaledkfreq, file = paste(tfname, '_', substr(names(scalefactors[j]), 1, nchar(names(scalefactors[j]))-27), '.rds', sep = ''))
       all_scalefactors[[j]] = scaledkfreq 
       }
      return(all_scalefactors)
}

##Take sum of each position's scaledkmerfreq (1/scale_factor * kmer frequency) and use to create 
##data.table of inputs for each mask for each TF motif
scaledkmerfreq.sum.plus <- function(input) {
    output = data.table(matrix(nrow = length(input[[1]])))
    for (i in 1:length(input)) {
        for (p in 1:length(input[[i]])) {
            output[p, names(input[i]) := sum(input[[i]][[p]][,4])]
        }
    }
    output[,1] <- NULL
    return(output)
}

scaledkmerfreq.sum.minus <- function(input) {
    output = data.table(matrix(nrow = length(input[[1]])))
    for (i in 1:length(input)) {
        for (p in 1:length(input[[i]])) {
            output[p, names(input[i]) := sum(input[[i]][[p]][,5])]
        }
    }
    output[,1] <- NULL
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
  output[, position := paste(rep('X',20,sep = ''), collapse = '')]
  
  for (i in 1:nrow(output)) {
    output[i,]$position = paste(substr(output[i, position],1,i-1),'N',substr(output[i, position],i+1,nchar(output[i, position])), sep = '')
    output[i,]$position = paste(substr(output[i, position],1,(nchar(output[i, position])/2)), 'C', 
                                paste(substr(output[i, position],(nchar(output[i, position])/2+1),nchar(output[i, position]))), sep = '')
  }
  return(output)
}


#


