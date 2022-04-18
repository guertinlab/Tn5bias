library(data.table)
library(matrixStats)
library(ggseqlogo)
library(lattice)
library(bigWig)
#FIMO.to.Bed takes input in FIMO format (1-based start and end) and converts 
#it to BED format (0-based start, 1-based end)
FIMO.to.BED <- function(fimofile) {
    require(data.table)
    #load in FIMO file
    FIMOdf <- fread(fimofile, sep = '\t', header = TRUE)
    #shift start position by -1 to make it 0-based
    FIMOdf[,4] = FIMOdf[,4] - 1
    #make location column
    FIMOdf[,location := paste(motif_id, ':', start, '-', stop, sep = '')]
    #Rearrange columns for BED format
    beddf <- FIMOdf[, c(3,4,5,2,11,6)]
    colnames(beddf) <- c('chr', 'start', 'end', 'gene', 'location', 'strand')
    return(beddf)
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
                                   factor, group, 
                                   stringsAsFactors = FALSE)
    colnames(composite.lattice) = c('x', 'est', 'factor', 'group')
    composite.lattice$x = as.numeric(composite.lattice$x)
    unload.bigWig(bw.plus)
    unload.bigWig(bw.minus)
    return(composite.lattice)
}
#

##Plot.composites takes the composite.lattice object and creates a plot for 
#this data, while also allowing a mapping of the original motif width onto the
#plot by specifying a TFlen value (the default is 10). 
plot.composites <- function(dat, ylabel = '', pdf_name = 'PLEASE_SET_FILE_NAME',
                            xlabel = '', striplabel = TRUE, legend = TRUE,
                            motifline = FALSE, Motiflen = 10,
                            figwidth = 2.5, figheight=3,
                            indexlist = NULL, layoutgrid = NULL,
                    col.lines = c("#0000FF", "#FF0000", "#008000",
                                  rgb(0,0,0,1/2),  
                                  rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), 
                                  rgb(1/2,1/2,0,1/2)), 
                                  fill.poly = c(rgb(0,0,1,1/4), 
                                  rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4), 
                                  rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    require(lattice)
    pdf(paste(pdf_name, '.pdf', sep = ''), width= figwidth, height= figheight) 
    print(xyplot(est ~ x|factor, group = group, data = dat, strip = striplabel,
                 type = 'l', as.table = TRUE, ylim=c(0, 0.02),
                 scales=list(x=list(cex=0.8,relation = "free", axs ="i"), 
                             y =list(cex=0.8, relation="free", tick.number=3)),
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
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 lwd=2,
                 ylab = list(label = paste(ylabel), cex =0.8),
                 xlab = list(label = paste(xlabel), cex =0.8),
                 index.cond = indexlist,
                 layout = layoutgrid,
                 panel = function(x, y, ...) {
                     panel.xyplot(x, y, ...)
                     #panel.abline(h = 0, lty =1, lwd = 1.0, col = '#A9A9A932')
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

