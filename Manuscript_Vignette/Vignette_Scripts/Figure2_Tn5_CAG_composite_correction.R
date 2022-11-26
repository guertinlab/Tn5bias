library(bigWig)
library(lattice)
library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
####CAG peak direction
#Load in random 400,000 CAG instances 
CAG_BED <- fread(input = 'CAG_locations_rand400k.bed')
colnames(CAG_BED) <- c('chr', 'start', 'end', 'gene', 'location', 'strand')

#Plot the unscaled signal at the CAG instances 
CAG_unscaled_composite <- BED.query.bigWig(CAG_BED,
                                           bwPlus = '../Figure1/C1_gDNA_rep1_plus.bigWig',
                                           bwMinus = '../Figure1/C1_gDNA_rep1_minus.bigWig',
                                           upstream = 15, downstream = 15, group = 'Unscaled', factor = '')

#Plot the NCNN scaled signal at the CAG instances 
CAG_NCNN_composite <- BED.query.bigWig(CAG_BED,
                                       bwPlus =  'C1_gDNA_rep1_plus_NCNN.bigWig',
                                       bwMinus =  'C1_gDNA_rep1_minus_NCNN.bigWig',
                                       upstream = 15, downstream = 15, group = 'NCNN', factor = '')
CAG_NCNN_composite <- rbind(CAG_NCNN_composite, CAG_unscaled_composite)
CAG_NCNN_composite$factor <- c('NCNN')


#Plot the CXXXXNNN scaled signal at the CAG instances 
CAG_CXXXXNNN_composite <- BED.query.bigWig(CAG_BED,
                                          bwPlus =  'C1_gDNA_rep1_plus_CXXXXNNN.bigWig',
                                          bwMinus =  'C1_gDNA_rep1_minus_CXXXXNNN.bigWig',
                                          upstream = 15, downstream = 15, group = 'CXXXXNNN', factor = '')
CAG_CXXXXNNN_composite <- rbind(CAG_CXXXXNNN_composite, CAG_unscaled_composite)
CAG_CXXXXNNN_composite$factor <- c('CXXXXNNN')

#Plot the NNNXXXC scaled signal at the CAG instances 
CAG_NNNXXXC_composite <- BED.query.bigWig(CAG_BED,
                                               bwPlus =  'C1_gDNA_rep1_plus_NNNXXXC.bigWig',
                                               bwMinus =  'C1_gDNA_rep1_minus_NNNXXXC.bigWig',
                                               upstream = 15, downstream = 15, group = 'NNNXXXC', factor = '')
CAG_NNNXXXC_composite <- rbind(CAG_NNNXXXC_composite, CAG_unscaled_composite)
CAG_NNNXXXC_composite$factor <- c('NNNXXXC')
CAG_fig_composite <- rbind(CAG_NCNN_composite,
                           CAG_CXXXXNNN_composite, CAG_NNNXXXC_composite)
save(CAG_fig_composite, file = 'Figure2B_CAG_correction.Rdata')
#Plot random CAG composites with mask corrections
plot.composites(CAG_fig_composite, legend = FALSE, 
                pdf_name = 'Figure2B_CAG_direction_maskcompare',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from CAG Center',
                indexlist = list(c(3,2,1)),
                layoutgrid = c(3,1),
                figwidth = 8, figheight=5,
                motifline = FALSE, y_axis = TRUE,
                y_axis_range = seq(0,0.0185, 0.003), x_axis_range = -15:15, nYticks = 2,
                col.lines = c('#096000', '#5409C9', '#FF0000','#A1A3AB'),
                linethick = 3)
