library(bigWig)
library(lattice)
library(data.table)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
load('Figure2B_CAG_correction.Rdata')
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
