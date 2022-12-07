source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)
library(lattice)
library(grid)
###########################
################################################################################
#Plot DNase unscaled composite overlay
#Change DNase group to a given factor, then change the factor column to 'DNase' 
load('../Figure5/DNase_unscaled_compositelist.Rdata')
DNase_unscaled_composites_overlay = do.call(rbind, DNase_unscaled_compositelist)
DNase_unscaled_composites_overlay$group = DNase_unscaled_composites_overlay$factor
DNase_unscaled_composites_overlay$factor = 'DNase'

plot.composites(DNase_unscaled_composites_overlay, legend = FALSE, 
                pdf_name = 'Figure7A_DNase_unscaled_overlay',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                Y_ticks = FALSE,col.lines = c("#8080804D"),
                y_axis_range = seq(0,0.05,0.02), y_axis = TRUE, labsize = 2.0)

################################################################################
#Plot Tn5 unscaled composite overlay
#Change Tn5 group to a given factor, then change the factor column to 'Tn5' 
load('Tn5_unscaled_compositelist.Rdata')  
Tn5_unscaled_composites_overlay = do.call(rbind, Tn5_unscaled_compositelist)
Tn5_unscaled_composites_overlay$group = Tn5_unscaled_composites_overlay$factor
Tn5_unscaled_composites_overlay$factor = 'Tn5'

plot.composites(Tn5_unscaled_composites_overlay, legend = FALSE, 
                pdf_name = 'Figure7A_Tn5_unscaled_overlay',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                Y_ticks = FALSE,col.lines = c("#8080804D"),
                y_axis_range = seq(0,0.05,0.02), y_axis = TRUE, labsize = 2.0)  
################################################################################
#Plot Tn5 unscaled E2F1 and Mef2A composites
Tn5_unscaled_compare_composites_overlay = do.call(rbind, Tn5_unscaled_compositelist[c(8,25)])
Tn5_unscaled_compare_composites_overlay$group = Tn5_unscaled_compare_composites_overlay$factor
Tn5_unscaled_compare_composites_overlay$factor = 'Tn5'

plot.composites(Tn5_unscaled_compare_composites_overlay, legend = FALSE, 
                pdf_name = 'Figure7C_Tn5_E2F1_MEF2A_unscaled_overlay',
                ylabel = '',
                xlabel = '',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                Y_ticks = FALSE,col.lines = c("#415287", "#F84C1E"),
                y_axis_range = seq(0,0.04,0.02), y_axis = TRUE,
                striplabel = FALSE, labsize = 2.0, linethick = 3)  
