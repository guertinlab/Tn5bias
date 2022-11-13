library(ggplot2)
load('../data/Figure5A_DNase_redundant_var_imps.Rdata')
redundant_base_imps$position = seq(-9.5,9.5,1)

setwd('../output')

pdf(file = 'Figure5A_Redundant_DNasePN_PRE_variable_importances.pdf', width = 9, height = 12)
impbarplot <-ggplot(data=redundant_base_imps, aes(x=position, y=possum)) 
impbarplot +
  labs(x= 'Position Relative to Cut Site', y= 'Importance') +
  geom_bar(stat="identity", fill="black") +
  theme_classic() +
  theme(axis.title.x = element_text(size=40, colour = "black", face = 'bold'),
        axis.title.y = element_text(size=40, colour = "black", face = 'bold'),
        axis.text.x = element_text(size=38, colour = "black"),
        axis.text.y = element_text(size=38, colour = "black"),
        axis.line = element_line(colour = 'black', size = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(-9.5, -4.5, 0.5, 4.5, 9.5), labels = c(-9.5, -4.5, 0.5, 4.5, 9.5)) + 
  theme(axis.text.x=element_text(vjust=1), plot.margin = margin(1,1,1.5,1.2, "cm"), axis.ticks.length=unit(.75, "cm"), axis.ticks=element_line(size=2))  
dev.off()

