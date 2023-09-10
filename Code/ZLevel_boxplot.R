### this code plots the discrimination of the outcome lables by each significant LF ###

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

x=read.csv("EvsNE_X.csv", row.names=1)
y=read.csv("EvsNE_Y.csv", row.names = 1)
x_std <- scale(x, T, T)

##loading er result
er_result=readRDS("/EvsNE.output/final_delta_0.1_lambda_1.rds")
##getting zs
z=read.csv("zmatrix_EvsNE.csv",row.names=1)


annot_colors=c("Egg+"="#E15759","Egg-"="#4E79A7")

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
file_name=paste0("Z16_boxplot_EvsNE.pdf")
file_name
ggsave(filename=file_name,plot=resplot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)
