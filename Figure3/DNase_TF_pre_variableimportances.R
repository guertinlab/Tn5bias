library(ggplot2)
load('DNase_redundant_var_imps.Rdata')

pdf(file = 'Redundant_DNasePN_PRE_variable_importances.pdf', width = 15, height = 20)
impbarplot <-ggplot(data=DNase_redundant_var_imps, aes(x=position, y=possum)) 
impbarplot +
  labs(x= 'Position Relative to Cut Site', y= 'Importance') +
  geom_bar(stat="identity", fill="black") +
  theme_classic() +
  theme(axis.title.x = element_text(size=36, face="bold", colour = "black"),
        axis.title.y = element_text(size=36, face="bold", colour = "black"),
        axis.text.x = element_text(size=26, face="bold", colour = "black"),
        axis.text.y = element_text(size=26, face="bold", colour = "black"),
        axis.line = element_line(colour = 'black', size = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(seq(-9.5, 9.5, 5)), labels = c(seq(-9.5, 9.5, 5))) + 
  theme(axis.text.x=element_text(vjust=1), plot.margin = margin(1,1,1.5,1.2, "cm"))  
dev.off()

