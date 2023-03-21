library(lattice)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
################################################################################
load('Figure7E_overlay_list.Rdata')
#Plot overlays:
################################################################################
################################################################################
#Plot seqOutBias overlay
plot.composites(Figure7E_overlay_list[[1]], legend = FALSE, 
                pdf_name = 'Figure7E_Tn5_seqOutBias_overlay',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                hline_val = 0.006140407, Y_ticks = FALSE,col.lines = c("#8080804D"),
                hline = TRUE, hlinelwd = 10, hlinelty = 3,
                y_axis_range = seq(0,0.05,0.02), y_axis = TRUE, labsize = 2.0)  
#Plot rule ensemble overlay
plot.composites(Figure7E_overlay_list[[2]], legend = FALSE, 
                pdf_name = 'Figure7E_Tn5_RE_overlay',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                hline_val = 0.006140407, Y_ticks = FALSE,col.lines = c("#8080804D"),
                hline = TRUE, hlinelwd = 10, hlinelty = 3,
                y_axis_range = seq(0,0.05,0.02), y_axis = TRUE, labsize = 2.0)  
#Plot unscaled overlay
plot.composites(Figure7E_overlay_list[[3]], legend = FALSE, 
                pdf_name = 'Figure7E_Tn5_unscaled_overlay',
                ylabel = 'Insertion Frequency',
                xlabel = 'Distance from Motif Center',
                motifline = FALSE, Motiflen = 0, figwidth = 12, figheight = 8,
                x_axis_range = -100:100, X_axis_ticks = seq(-100,100,50), 
                hline_val = 0.006140407, Y_ticks = FALSE,col.lines = c("#8080804D"),
                hline = TRUE, hlinelwd = 10, hlinelty = 3,
                y_axis_range = seq(0,0.05,0.02), y_axis = TRUE, labsize = 2.0)  
