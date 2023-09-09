## this code makes box plots for the results obtained from LASSO_LepnoLep_part1.R ##

library(ggplot2)
library(matrixStats)
library(ggedit)
library(grid)
library(gtools)
library(readxl)

sig_text <- "P<0.01"
sig_grob = grid.text(sig_text, x=0.5,  y=0.95, gp=gpar(col="black", fontsize=8, fontface="bold"))

not_sig_text <- "n.s."
not_sig_grob = grid.text(not_sig_text, x=0.5,  y=0.95, gp=gpar(col="black", fontsize=8, fontface="bold"))

dF <- read.csv("AUC_class2_SchistoAg0.7_RF.csv")
Actual <- as.vector(dF[1:10, 2])
Permuted <- as.vector(dF[11:20, 2])
wilcox.test(Actual, Permuted)
##median(Actual$AUC)
##median(Permuted$AUC)
dF$Type <- as.factor(dF$Type)
p <- ggplot(dF, aes(x=Type, y=AUC, colour = Type)) + geom_boxplot(width = 0.15) 
p + theme_classic()  + annotation_custom(sig_grob) + theme(legend.position="None")+ ylim(0.4, 1.0) + xlab(" ") + ylab("AUC_RF")

ggsave("AUC_LASSO.class2_SchistoAg0.7_RF.pdf")