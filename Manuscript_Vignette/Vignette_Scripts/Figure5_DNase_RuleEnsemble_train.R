library(data.table)
library(pre)
library(ggplot2)
source('https://raw.githubusercontent.com/guertinlab/Tn5bias/master/Manuscript_Vignette/Vignette_Scripts/Tn5_Bias_Functions.R')
options(scipen = 100)

#These are the training and test TF motif sets
train = c( 'ASCL1', 'ATF3', 'CLOCK', 'CTCF', 'DLX2', 'EBF1',
           'EGR1', 'ELF1', 'FOS', 'FOXA1', 'FOXK2', 'GATA2',
           'GLIS1', 'IRF1', 'JUN', 'LEF1', 'MEIS2', 'MLX',
           'MYC', 'PPARG','RUNX1', 'SPIB', 'SRY', 'STAT1', 'TGIF1')

test = c('AR', 'CEBPB', 'E2F1', 'ESR1', 'FERD3L', 'HOXC12',
         'HSF1', 'MAX', 'MEF2A', 'MGA', 'NFATC3', 'NR2F2',
         'POU3F1', 'REST', 'Six3', 'SP1', 'USF1', 'TEAD1')

load('DNase_TF_plus_pre_input.Rdata')
load('DNase_TF_minus_pre_input.Rdata')

#Subset test data


DNase_pre_data <- rbind(Plus_pred_input, minus_pred_input)

#Remove X, factor, group, correction columns
DNase_pre_data = data.frame(DNase_pre_data[,-c(1)])
################################################################################
#Train rule ensemble model using DNase input
set.seed(42)
DNase_pre_model <- pre(unscaled ~  ., data = DNase_pre_data, type = 'both', sampfrac = 0.5, 
                       verbose = TRUE, ntrees = 3500,
                       use.grad = TRUE, nfolds = 10, winsfrac = 0, maxdepth = 3L,
                       learnrate = 0.1, intercept = FALSE)
save(DNase_pre_model, file = 'DNase_PRE_model.Rdata')
################################################################################
#Write table for rule ensemble model coefficients/rules
print(DNase_pre_model, penalty.par.val = DNase_pre_model$glmnet.fit$lambda.min)
DNase_pre_model_coef <- coef(DNase_pre_model, penalty.par.val = DNase_pre_model$glmnet.fit$lambda.min)
DNase_pre_model_coef = DNase_pre_model_coef[which(DNase_pre_model_coef$coefficient != 0),]
write.table(DNase_pre_model_coef[,2:3], file = 'DNase_pre_coeff.txt', quote = FALSE, sep = ',',
            row.names = FALSE)
#Plot variable importances for PRE model:
imps = importance(DNase_pre_model, penalty.par.val = DNase_pre_model$glmnet.fit$lambda.min,
                  abbreviate = FALSE)
base_imps = imps$varimps
base_imps$varname = as.factor(base_imps$varname)
redundant_base_imps = redundantpos(base_imps, labelnum = 1, impnum = 2)
redundant_base_imps = redundant_base_imps[9:28,]
redundant_base_imps$position = seq(-9.5,9.5,1)
#Plot redundant_var_imps
pdf(file = 'Figure5A_Redundant_DNase_PRE_variable_importances.pdf', width = 9, height = 12)
impbarplot <-ggplot(data=redundant_base_imps, aes(x=position, y=possum)) 
impbarplot +
  labs(x= 'Position Relative to Cut Site', y= 'Importance') +
  geom_bar(stat="identity", fill="black") +
  theme_classic() +
  theme(axis.title.x = element_text(size=40, colour = "black"),
        axis.title.y = element_text(size=40, colour = "black"),
        axis.text.x = element_text(size=38, colour = "black"),
        axis.text.y = element_text(size=38, colour = "black"),
        axis.line = element_line(colour = 'black', size = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(-9.5, -4.5, 0.5, 4.5, 9.5), labels = c(-9.5, -4.5, 0.5, 4.5, 9.5)) + 
  scale_y_continuous(breaks = seq(0,0.6,0.1), labels = seq(0,0.6,0.1)) + 
  theme(axis.text.x=element_text(vjust=1), plot.margin = margin(1,1,1.5,1.2, "cm"),
        axis.ticks.length=unit(.75, "cm"), axis.ticks=element_line(size=2))  
dev.off()

