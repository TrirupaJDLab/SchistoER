### this code makes box plots for the results obtained from LASSO_LepnoLep_part1.R ###

##loading libraries 
library(ggplot2)
library(matrixStats)
library(grid)
library(gtools)
library(readxl)
library(tidyr)
library(dplyr)

## concatenating the real and permuted datafraes for each
RF_real <- read.csv("/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/aucRF_allpat_AbScAg_MFI50_l0.5.csv", header = TRUE)
RF_permute <- read.csv("/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/aucRFpermute_noLepPat_AbScAg_MFI50_l0.5.csv", header = TRUE)

#SVM_real= read.csv("/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/aucSVM_allpat_AbScAg_MFI50_l0.5.csv", header = TRUE)
#SVM_permute= read.csv("/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/aucSVMpermute_noLepPat_AbScAg_MFI50_l0.5.csv", header = TRUE)

# Assuming both CSV files have the same column names and order
RF_csv = rbind(RF_real, RF_permute)
#SVM_csv= rbind(SVM_real, SVM_permute)


#RF_csv=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/CSVfiles_v2/maxExp.sub.var25.cov1.RF_auc_0.5.ntree15_stitched.csv")
#SVM_csv=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/CSVfiles_v2/maxExp.sub.var25.cov1.SVM_auc_0.5.cost20_stitched.csv")

##changing column names
RF_df=as.data.frame(RF_csv[,-1]) ###removing the first column that has repCV number
#SVM_df=as.data.frame(SVM_csv[,-1]) ## also removing the 21st element since it was the column header while concatenating actual +permuted files

##now we need to label the Actual and Permuted data points as a second column
Type_actual <-rep("Actual",times=20)
Type_permuted =rep("Permuted",times=20)
RF_df$Type=c(Type_actual,Type_permuted)
SVM_df$Type=c(Type_actual,Type_permuted)
colnames(RF_df)[1]<- "AUC"
colnames(SVM_df)[1]<- "AUC"
RF_df$AUC=as.numeric(RF_df$AUC)
SVM_df$AUC=as.numeric(SVM_df$AUC)
# setwd("/Users/sar210/Box/MSD_data_ABMR/")
#setwd("/ix/djishnu/Trirupa/ABomics.Prj/Lasso/miRNA/")


sig_text <- "P<0.01"
sig_grob = grid.text(sig_text, x=0.5,  y=0.95, gp=gpar(col="black", fontsize=8, fontface="bold"))

not_sig_text <- "n.s."
not_sig_grob = grid.text(not_sig_text, x=0.5,  y=0.95, gp=gpar(col="black", fontsize=8, fontface="bold"))

##AUCs for RF
#significance test for RF_df 
Actual <- as.vector(RF_df[1:20, 1])
Permuted <- as.vector(RF_df[21:40, 1])
wilcox.test(Actual, Permuted)
median(Actual)
median(Permuted)
##plotting box plot for RF_df
RF_df$Type <- as.factor(RF_df$Type)
p <- ggplot(RF_df, aes(x=Type, y=AUC, colour = Type)) + geom_boxplot(width = 0.15) 
p + theme_classic()  + annotation_custom(sig_grob) + theme(legend.position="None")+ ylim(0.4, 0.9) + xlab(" ") + ylab("AUC_RF")

ggsave("/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/boxplotAUC_RF_allpat_AbScAg_MFI50_30May23_l0.5.pdf",device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)
