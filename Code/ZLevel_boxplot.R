### this code plots the discrimination of the outcome lables by each significant LF ###
#rm(list = ls())
#cat("\014")

library(EssReg)
library(doParallel)
library(dplyr)
library(pROC)
library(ROCR)
library(SLIDE)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(reshape2)
registerDoParallel(detectCores())

x=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/EvsNE_AbScAg_30May23_MFI50_X.csv", row.names=1)
y=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/EvsNE_Ab_ScAg_30May23_Y.csv", row.names = 1)
x_std <- scale(x, T, T)

##loading er result
er_result=readRDS("final_delta_0.1_lambda_1.rds")
##getting zs
z=read.csv("zmatrix_EvsNE.MFI50.d0.1_lam1.auc.csv",row.names=1)

sig_idx <- c(10)
sig_z <- as.data.frame(z[ ,10])
print(colnames(sig_z))

##renaming Ys
y$Y[y$Y == 1] <- "EndoA"
y$Y[y$Y == 0] <- "EndoB"

## for significance calculation
signif_labels=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

df1 <- data.frame()
df1=cbind(sig_z,y[,1])
colnames(df1)[ncol(df1)] ="Status"
colnames(df1)[1] ="Z10"




annot_colors=c("Egg+"="#E15759","Egg-"="#4E79A7")
#annot_colors=c("SEA+Egg+"="#E15759","SEA+Egg-"="#F28E2B")
#annot_colors=c("SEA+Egg+"="#E15759","EndoB"="#FFBF00")
##renaming Ys
y$Y[y$Y == 1] <- "Egg+"
y$Y[y$Y == 0] <- "Egg-"
sig_z <- as.data.frame(z[ ,16])
print(colnames(sig_z))
i="Z16"
colnames(sig_z)[1] ="Z16"  ##Z12=protective calu LF
value <- sig_z[ ,i]
df1 <- as.data.frame(value)
df1["Status"] <- y[, 1]

signif_labels=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

df1_zplot=ggplot2::ggplot(data=df1,aes(x=Status,y=value))+
  ggplot2::geom_boxplot(width=0.40, aes(fill=as.factor(Status))) + 
  scale_fill_manual(values=annot_colors)+
  ggtitle("Glycodominant LF") +theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), panel.spacing = unit(0.5, "lines"),
                    panel.border = element_rect(fill = NA, color = "black",linetype="dashed"),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_compare_means(symnum.args = signif_labels,comparison=list(c("Egg+","Egg-")),
                    aes(label = paste0("p = ", after_stat(p.format))),label.x.npc="middle",label.y.npc="top")

df1_zplot                                   
#resplot=df1_zplot+ylim(min(df1$value),max(df1$value))
resplot=df1_zplot+ylim(min(df1$value),4)
resplot=df1_zplot
resplot
file_name=paste0("/ix/djishnu/Trirupa/Schisto_Proj/ER_run/Egg_NoEgg_revisedLabels/Plots/version2_30May23/Z16_boxplot_EvsNE_30May23.pdf")
file_name
ggsave(filename=file_name,plot=resplot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)
