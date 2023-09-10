#### This code gets AUCs of features picked up by ER selected LFs in both discovery and validation cohorts ### 

library(EssReg)
library(doParallel)
library(dplyr)
library(pROC)
library(ROCR)
library(SLIDE)
library(matrixStats)
library(ggplot2)
library(glmnet)
registerDoParallel(detectCores())

x=read.csv("active_endoB_X.csv", row.names=1)
y=read.csv("active_endoB_Y.csv", row.names = 1)
x_std <- scale(x, T, T)

##loading er result
er_result=readRDS("./active_endoB.output/final_delta_0.1_lambda_1.rds")
##getting zs
z <- predZ(x_std, er_result)
colnames(z)=paste0("Z",seq(1,ncol(z)))
z=read.csv("./active_endoB.output/zmatrix_active_endoB.csv",row.names=1)

##getting features within marginals
readRDS(file="./active_endoB.output/hierER.active_endoB.spec0.1.rds")
A_abs <-  apply(er_result$A,c(1,2),abs)
A_k <- as.data.frame(A_abs) %>% select(12) ##insert LF number
var_name <- rownames(subset(A_k,A_k > 0.08)) 
var_name

result_df=data.frame(name=NULL,pval=NULL, aucXiY=NULL,corXiY=NULL,A_value=NULL)

for (i in var_name){
  res=wilcox.test(x[,i]~as.matrix(y));
  auc_score=(pROC:::roc(as.factor(y$Y),x[,i], direction="<"))$auc
  result_df<-rbind(result_df,data.frame(name=i,pval=res$p.value,corXiY=cor(x[,i],y),aucXiY=auc_score , A_value=A_k[i,]))
  }
sorted_result_df <- result_df[order(result_df$A_value,decreasing=TRUE),]
colnames(sorted_result_df)<- c("name","pval","corrXiY","aucXiY","A_value")
rownames(sorted_result_df)=sorted_result_df$name
final_df=sorted_result_df[,2:ncol(sorted_result_df)]


### checking if aucs from the Validation cohort are also tracking with Y similarly
VC_data=read.csv("XY_VC_imputed.csv")
Y=VC_data$Y

##ER selected features
Z1=c("IgG2.Gst.SM25", "IgG2.Calumenin", "IgG3.SEA", "IgG3.CD63", "IgG4.SEA", "IgG4.CD63", "IgA1.MEG", "FcREpsilonI.SEA", "FcR2A.Calumenin")
Z2=c("IgG.Calumenin", "IgG3.Gst.SM25", "IgA1.SEA", "IgA1.Gst.SM25", "IgA1.MEG", "IgA1.CD63", "IgA1.Calumenin")
Z3=c("IgG.Calumenin", "IgG1.Gst.SM25", "IgG1.Calumenin", "FcR1.CD63", "FcR1.Calumenin")
Z4=c("IgG.Calumenin", "IgG1.CD63", "IgG1.Calumenin", "IgG2.Calumenin", "FcR2A.Calumenin", "FcRIIB.Calumenin", "FcRIIIA.Calumenin", "FcR3B.Calumenin")
Z5=c("IgM.MEG", "FcR2A.SEA", "RCA.MEG", "SNA.Gst.SM25", "SNA.MEG", "SNA.CD63")
allZfeats=c(Z1,Z2,Z3,Z4,Z5)
ERfeats=unique(allZfeats)

auc_df=data.frame(name=NULL,aucXiY=NULL)
for (i in ERfeats){
  auc_score=(pROC:::roc(as.factor(VC_data$Y),VC_data[,i], direction="<"))$auc
  auc_df<-rbind(auc_df,data.frame(name=i,aucXiY=auc_score))
}
auc_df

write.csv(auc_df,"ERselectedfeat_directionsInVC.csv")
