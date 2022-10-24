library(ggplot2)
load('../data/Figure6C_Tn5_RE_redundant_var_imps.Rdata')
redundant_var_imps$position = seq(-19,19,1) 

setwd('../output')

pdf(file = 'Figure6C_Redundant_Tn5_RE_variable_importances.pdf', width = 12, height = 12)
impbarplot <-ggplot(data=redundant_var_imps, aes(x=position, y=possum)) 
impbarplot +
  labs(x= 'Position Relative to Central Base', y= 'Importance') +
  geom_bar(stat="identity", fill="black") +
  theme_classic() +
  theme(axis.title.x = element_text(size=40, colour = "black", face = 'bold'),
        axis.title.y = element_text(size=40, colour = "black", face = 'bold'),
        axis.text.x = element_text(size=38, colour = "black"),
        axis.text.y = element_text(size=38, colour = "black"),
        axis.line = element_line(colour = 'black', size = 1),
        plot.title = element_text(hjust = 0.5)) +
  #scale_x_continuous(breaks = c(-9.5, -4.5, 0.5, 4.5, 9.5), labels = c(-9.5, -4.5, 0.5, 4.5, 9.5)) + 
  theme(axis.text.x=element_text(vjust=1), plot.margin = margin(1,1,1.5,1.2, "cm"), axis.ticks.length=unit(.75, "cm"), axis.ticks=element_line(size=2)) +
  scale_x_continuous(breaks = c(-15,0,15))
dev.off()

