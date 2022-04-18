library(lattice)
source('composite_functions.R')


load('Figure2B_Tn5_CAG_peak_direction_composite.Rdata')

plot.composites(Tn5_CAG_peak_direction_composite, legend = FALSE, 
                pdf_name = 'Figure2B_CAG_direction_maskcompare',
                ylabel = 'Cut Frequency',
                xlabel = 'Distance from Motif Center',
                indexlist = list(c(3,2,1)),
                layoutgrid = c(3,1),
                figwidth = 9, figheight=5,
                motifline = FALSE)
