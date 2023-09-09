## this code gets the signficant LFs from the output of ER pipeline and the features within it. ##

library(EssReg)
library(doParallel)
library(dplyr)
library(pROC)
library(ROCR)
library(SLIDE)
library(matrixStats)
library(ggplot2)
registerDoParallel(detectCores())

x=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/unsupV2/class2_AbScAg.MFI50_30May23_V2_X.csv", row.names=1)
y=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/unsupV2/class2_AbScAg.MFI50_30May23_V2_Y.csv", row.names = 1)
x_std <- scale(x, T, T)

##loading er result
er_result=readRDS("final_delta_0.1_lambda_1.rds")
##getting zs
z <- predZ(x_std, er_result)
colnames(z)=paste0("Z",seq(1,ncol(z)))

#write.csv(z,file="zmatrix_c2.MFI50.d0.05_lam1.auc.csv",row.names=TRUE)
z=read.csv("zmatrix_c2.MFI50.d0.1_lam1.auc.csv",row.names=1)

## running SLIDE (without interaction terms) ##significant latent factors
SLIDE_marginal_res =SLIDE(z,y,niter=1000,do_interacts = F, spec=0.2, f_size = 39,fdr = 0.1)
SLIDE_marginal_res
readRDS(file="hierER.c2.MFI50_SchistoAg.d0.1_lam1_spec0.2.rds")

### getting the features at an A threshold, from each LF
A_abs <-  apply(er_result$A,c(1,2),abs)
A_k <- as.data.frame(A_abs) %>% select(10)
var_name <- rownames(subset(A_k,A_k > 0.08)) 
var_name