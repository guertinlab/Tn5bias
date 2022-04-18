library(ggplot2)
################################################################################
################################################################################
################################################################################


load('Figure3B_DNase_PRE_NNNCNN_log_fold_change.Rdata')

Maskvarplot <- function(){
  ggplot(DNase_PRE_NNNCNN_log_fold_change, aes(x=factor(Factor), y=Difference, fill = Treatment, color = Treatment)) +
    xlab('') + ylab(expression(log[2]*('Fold Change'))) +
    geom_violin() +
    scale_fill_manual(values=c('#808080', '#FF0000', '#0000FF')) +
    scale_color_manual(values=c('#808080', '#FF0000', '#0000FF')) +
    #scale_x_discrete(position = "top") +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #geom_point(size = 1, shape = 20) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75, size = 16, face = 'bold', color = 'black'), 
          axis.title.y=element_text(size=16, color = 'black', face = 'bold'),
          axis.title.x = element_text(vjust = -2),
          axis.text.y = element_text(size = 8, color = 'black', face = 'bold'),
          axis.line = element_line(colour = 'black', size = 1),
  legend.title = element_text(size=0),
  legend.text = element_text(size=16, face = 'bold'))
}

pdf("Figure3B_DNase_PRE_NNNCNN_testset_log2.pdf", width=12, height=7)
Maskvarplot()
dev.off()

