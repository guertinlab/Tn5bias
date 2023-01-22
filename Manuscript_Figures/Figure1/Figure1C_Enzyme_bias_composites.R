source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
################################################################################
load('../data/Figure1C_Unscaled_enzyme_composites.Rdata')
###Plot these all together
enzymebiasplot = plot.composites(unscaled_composite_agg, x_axis_range=-15:15,
                                 legend = FALSE, 
                                 ylabel = 'Cut Frequency',
                                 xlabel = 'Distance from Motif Center',
                                 figwidth = 4, figheight = 12,
                                 motifline = FALSE,
                                 Motiflen = mlen,
                                 striplabel = TRUE,
                                 pdf_name = 'unscaled_enzymebias_composites',
                                 indexlist = list(c(5,3,2,1,4)),
                                 layoutgrid = c(1,5),
                                 y_axis_range = 0:1, y_axis = TRUE, Y_axis_ticks = seq(0,1,0.5),
)
