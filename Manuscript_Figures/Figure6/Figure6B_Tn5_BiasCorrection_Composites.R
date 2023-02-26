library(lattice)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
################################################################################
load('Figure6B_Tn5_composites.Rdata')
load('supplemental_Figure6B_Tn5_composites.Rdata')
#Plot composites
plot.composites(Figure6B_list[[1]], legend = TRUE, 
                pdf_name = 'Figure6B_Tn5_RuleEnsemble_seqOutBias_comparison',
                figwidth = 8, figheight = 4,
                ylabel = '',
                xlabel = '',
                motifline = FALSE, Motiflen = Figure6B_list[[2]], layoutgrid = c(3,1),
                x_axis_range = -20:20, X_axis_ticks = seq(-20,20,10),
                Y_axis_ticks = seq(0,0.05,0.01), nYaxisdigits = 2,
                hline_val = 0.006140407, y_axis_range = seq(0,0.02,0.01),
                y_axis = FALSE, hline = TRUE, Y_ticks = FALSE, labsize = 0.85)

#Plot supplemental composites
plot.composites(supp_Figure6B_list[[1]], legend = TRUE, 
                pdf_name = 'Supplemental_Tn5_RuleEnsemble_seqOutBias_comparison',
                figwidth = 12, figheight = 10,
                ylabel = '',
                xlabel = '',
                motifline = TRUE, Motiflen = supp_Figure6B_list[[2]], layoutgrid = c(5,3),
                x_axis_range = -20:20, X_axis_ticks = seq(-20,20,10),
                Y_ticks = FALSE,
                hline_val = 0.006140407, y_axis_range = seq(0,0.02,0.01),
                y_axis = FALSE, hline = TRUE, labsize = 0.85)
