## this code is to do LASSO regression using ER selected features to show with Antigen-based models are better than SEA in Schisto prediction
#rm(list = ls())
#cat("\014")

library(e1071)
library(caret)
library(glmnet)
library(randomForest)
library(tree)
library(gbm)
library(matrixStats)
library(readxl)
library(cvAUC)
library(pROC)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

## loading datasets
Xdf=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/unsupV2/class2_AbScAg.MFI50_30May23_V2_X.csv", row.names=1)
Ydf=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/unsupV2/class2_AbScAg.MFI50_30May23_V2_Y.csv", row.names=1)

#ag_lis=c("SEA","CD63","Calumenin")
## subsetting only the selected Calumenin features for c1 vs c2
# Subset columns containing "Calumenin"
feats_Xdf <- Xdf[, grepl("CD63", colnames(Xdf))]


#### LASSO regression 
##perform k-fold cross-validation to find optimal lambda value
######## Splitting data into test and train ######

num_repCV=20
rep_aucs_CD63=c()
for (replicateCV in 1:num_repCV){
  
  k=5
  print(paste0("rep:",replicateCV))
  folds <- createFolds(y = Ydf$Y, k=k, list = FALSE, returnTrain = FALSE)
  # at each fold, we make sure to start with the true y label
  Y=Ydf$Y
  ori_myData <- cbind(Y,feats_Xdf, folds)
  fold_aucs=c()
  #for each fold of the data
  for (NoF in 1:k){
  fold = which(folds == NoF)
  #print(dim(ori_myData))
  myData <- ori_myData
  
  train <- myData[-fold, ]
  test <- myData[fold, ]
  
  # delete the fold column
  train <- train[, -ncol(train)]
  test <- test[, -ncol(test)]

  xtrain = as.matrix(train[, -1])
  ytrain=train$Y
  
  xtest=as.matrix(test[,-1])
  ytest=test$Y
  
  #### move to next repCV if only one class of ytest is present ####
  if (length(unique(ytest)) != 2){next}
  ######### train model #################
  cv_model <- cv.glmnet(x=xtrain, y=ytrain, alpha = 1, nfold=10)

  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  best_lambda
  tuned_lambda <-best_lambda*0.8    ##multiplication factor=0.7 (c2vc3), 0.8 for EvsNE
  #produce plot of test MSE by lambda value
  #plot(cv_model)

  #find coefficients of best model
  best_model <- glmnet(x=xtrain, y=ytrain, alpha = 1, lambda = tuned_lambda)
  coef(best_model)

  ########### prediction and validation on test data #################
  y_pred=predict(best_model, s=tuned_lambda, newx=xtest) ##this is the log-odds (of it being either 0 or 1) ratio 
  auc_pred=auc(ytest,y_pred)    
  fold_aucs=append(fold_aucs,auc_pred)
  }
  print(fold_aucs)
  median_fold_aucs=median(fold_aucs)
  rep_aucs_CD63=append(rep_aucs_CD63,median_fold_aucs)
  
}
median(rep_aucs_CD63)

### doing the same for SEA
feats_Xdf <- Xdf[, grepl("SEA", colnames(Xdf))]

num_repCV=20
rep_aucs_SEA=c()

for (replicateCV in 1:num_repCV){
  k=5
  print(paste0("rep:",replicateCV))
  folds <- createFolds(y = Ydf$Y, k=k, list = FALSE, returnTrain = FALSE)
  # at each fold, we make sure to start with the true y label
  Y=Ydf$Y
  ori_myData <- cbind(Y,feats_Xdf, folds)
  fold_aucs=c()
  #for each fold of the data
  for (NoF in 1:k){
    fold = which(folds == NoF)
    #print(dim(ori_myData))
    myData <- ori_myData
    
    train <- myData[-fold, ]
    test <- myData[fold, ]
    
    # delete the fold column
    train <- train[, -ncol(train)]
    test <- test[, -ncol(test)]
    
    xtrain = as.matrix(train[, -1])
    ytrain=train$Y
    
    xtest=as.matrix(test[,-1])
    ytest=test$Y
    
    #### move to next repCV if only one class of ytest is present ####
    if (length(unique(ytest)) != 2){next}
    ######### train model #################
    cv_model <- cv.glmnet(x=xtrain, y=ytrain, alpha = 1, nfold=10)
    
    #find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min
    best_lambda
    tuned_lambda=best_lambda*0.8 ## 0.7 for EvsNE # 0.75 for endo2b vs c3 # 0.9 for endoA vs c3
    #produce plot of test MSE by lambda value
    #plot(cv_model)
    
    #find coefficients of best model
    best_model <- glmnet(x=xtrain, y=ytrain, alpha = 1, lambda = tuned_lambda)
    coef(best_model)
    
    ########### prediction and validation on test data #################
    y_pred=predict(best_model, s=tuned_lambda, newx=xtest) ##this is the log-odds (of it being either 0 or 1) ratio 
    auc_pred=auc(ytest,y_pred)    
    fold_aucs=append(fold_aucs,auc_pred)
  }
  print(fold_aucs)
  median_fold_aucs=median(fold_aucs)
  rep_aucs_SEA=append(rep_aucs_SEA,median_fold_aucs)
  
}
median(rep_aucs_SEA)

# Subset columns containing "Calumenin"
feats_Xdf <- Xdf[, grepl("Calumenin", colnames(Xdf))]

#### LASSO regression 
##perform k-fold cross-validation to find optimal lambda value
######## Splitting data into test and train ######

num_repCV=20
rep_aucs_Cal=c()
for (replicateCV in 1:num_repCV){
  
  k=5
  print(paste0("rep:",replicateCV))
  folds <- createFolds(y = Ydf$Y, k=k, list = FALSE, returnTrain = FALSE)
  # at each fold, we make sure to start with the true y label
  Y=Ydf$Y
  ori_myData <- cbind(Y,feats_Xdf, folds)
  fold_aucs=c()
  #for each fold of the data
  for (NoF in 1:k){
    fold = which(folds == NoF)
    #print(dim(ori_myData))
    myData <- ori_myData
    
    train <- myData[-fold, ]
    test <- myData[fold, ]
    
    # delete the fold column
    train <- train[, -ncol(train)]
    test <- test[, -ncol(test)]
    
    xtrain = as.matrix(train[, -1])
    ytrain=train$Y
    
    xtest=as.matrix(test[,-1])
    ytest=test$Y
    
    #### move to next repCV if only one class of ytest is present ####
    if (length(unique(ytest)) != 2){next}
    ######### train model #################
    cv_model <- cv.glmnet(x=xtrain, y=ytrain, alpha = 1, nfold=10)
    
    #find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min
    best_lambda
    tuned_lambda <-best_lambda*0.8  ## 0.8 for EvsNE
    #produce plot of test MSE by lambda value
    #plot(cv_model)
    
    #find coefficients of best model
    best_model <- glmnet(x=xtrain, y=ytrain, alpha = 1, lambda = tuned_lambda)
    coef(best_model)
    
    ########### prediction and validation on test data #################
    y_pred=predict(best_model, s=tuned_lambda, newx=xtest) ##this is the log-odds (of it being either 0 or 1) ratio 
    auc_pred=auc(ytest,y_pred)    
    fold_aucs=append(fold_aucs,auc_pred)
  }
  print(fold_aucs)
  median_fold_aucs=median(fold_aucs)
  rep_aucs_Cal=append(rep_aucs_Cal,median_fold_aucs)
  
}
median(rep_aucs_Cal)

# Subset columns containing "Calumenin"
feats_Xdf <- Xdf[, grepl("MEG", colnames(Xdf))]

#### LASSO regression 
##perform k-fold cross-validation to find optimal lambda value
######## Splitting data into test and train ######

num_repCV=20
rep_aucs_MEG=c()
for (replicateCV in 1:num_repCV){
  
  k=5
  print(paste0("rep:",replicateCV))
  folds <- createFolds(y = Ydf$Y, k=k, list = FALSE, returnTrain = FALSE)
  # at each fold, we make sure to start with the true y label
  Y=Ydf$Y
  ori_myData <- cbind(Y,feats_Xdf, folds)
  fold_aucs=c()
  #for each fold of the data
  for (NoF in 1:k){
    fold = which(folds == NoF)
    #print(dim(ori_myData))
    myData <- ori_myData
    
    train <- myData[-fold, ]
    test <- myData[fold, ]
    
    # delete the fold column
    train <- train[, -ncol(train)]
    test <- test[, -ncol(test)]
    
    xtrain = as.matrix(train[, -1])
    ytrain=train$Y
    
    xtest=as.matrix(test[,-1])
    ytest=test$Y
    
    #### move to next repCV if only one class of ytest is present ####
    if (length(unique(ytest)) != 2){next}
    ######### train model #################
    cv_model <- cv.glmnet(x=xtrain, y=ytrain, alpha = 1, nfold=10)
    
    #find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min
    best_lambda
    tuned_lambda <-best_lambda*1
    #produce plot of test MSE by lambda value
    #plot(cv_model)
    
    #find coefficients of best model
    best_model <- glmnet(x=xtrain, y=ytrain, alpha = 1, lambda = tuned_lambda)
    coef(best_model)
    
    ########### prediction and validation on test data #################
    y_pred=predict(best_model, s=tuned_lambda, newx=xtest) ##this is the log-odds (of it being either 0 or 1) ratio 
    auc_pred=auc(ytest,y_pred)    
    fold_aucs=append(fold_aucs,auc_pred)
  }
  print(fold_aucs)
  median_fold_aucs=median(fold_aucs)
  rep_aucs_MEG=append(rep_aucs_MEG,median_fold_aucs)
  
}
median(rep_aucs_MEG)

# Subset columns containing "Calumenin"
feats_Xdf <- Xdf[, grepl("SM25", colnames(Xdf))]

#### LASSO regression 
##perform k-fold cross-validation to find optimal lambda value
######## Splitting data into test and train ######

num_repCV=20
rep_aucs_SM25=c()
for (replicateCV in 1:num_repCV){
  
  k=5
  print(paste0("rep:",replicateCV))
  folds <- createFolds(y = Ydf$Y, k=k, list = FALSE, returnTrain = FALSE)
  # at each fold, we make sure to start with the true y label
  Y=Ydf$Y
  ori_myData <- cbind(Y,feats_Xdf, folds)
  fold_aucs=c()
  #for each fold of the data
  for (NoF in 1:k){
    fold = which(folds == NoF)
    #print(dim(ori_myData))
    myData <- ori_myData
    
    train <- myData[-fold, ]
    test <- myData[fold, ]
    
    # delete the fold column
    train <- train[, -ncol(train)]
    test <- test[, -ncol(test)]
    
    xtrain = as.matrix(train[, -1])
    ytrain=train$Y
    
    xtest=as.matrix(test[,-1])
    ytest=test$Y
    
    #### move to next repCV if only one class of ytest is present ####
    if (length(unique(ytest)) != 2){next}
    ######### train model #################
    cv_model <- cv.glmnet(x=xtrain, y=ytrain, alpha = 1, nfold=10)
    
    #find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min
    best_lambda
    tuned_lambda <-best_lambda*1 ## 0.7 for c2vsc3
    #produce plot of test MSE by lambda value
    #plot(cv_model)
    
    #find coefficients of best model
    best_model <- glmnet(x=xtrain, y=ytrain, alpha = 1, lambda = tuned_lambda)
    coef(best_model)
    
    ########### prediction and validation on test data #################
    y_pred=predict(best_model, s=tuned_lambda, newx=xtest) ##this is the log-odds (of it being either 0 or 1) ratio 
    auc_pred=auc(ytest,y_pred)    
    fold_aucs=append(fold_aucs,auc_pred)
  }
  print(fold_aucs)
  median_fold_aucs=median(fold_aucs)
  rep_aucs_SM25=append(rep_aucs_SM25,median_fold_aucs)
  
}
median(rep_aucs_SM25)


median(rep_aucs_SEA)
median(rep_aucs_CD63)
median(rep_aucs_Cal)
## for comparing between Antigen1 and SEA
Ag_aucdf=as.data.frame(rep_aucs_SEA,col.names=c("SEA"))
Ag_aucdf['SM25']=rep_aucs_SM25
Ag_aucdf['CD63']=rep_aucs_CD63
Ag_aucdf["Calumenin"]=rep_aucs_Cal
colnames(Ag_aucdf)[1]="SEA"
melted_Ag_aucdf = melt(Ag_aucdf)

signif_labels=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
Ag_plot=ggplot2::ggplot(data=melted_Ag_aucdf,aes(x=variable,y=value))+ggtitle("EndoA vs EndoB")+
  ggplot2::geom_boxplot(width=0.40, aes(fill=as.factor(variable)))+xlab("Antigen")+ylab("AUC")+
  
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.spacing = unit(0.5, "lines"),
        panel.border = element_rect(fill = NA, color = "black",linetype="dashed"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  stat_compare_means(comparisons=list(c("SEA","SM25")),symnum.args = signif_labels,aes(label = paste0("p = ", after_stat(p.format))),
  #stat_compare_means(comparisons=list(c("SEA","CD63"),c("SEA","Calumenin"),c("CD63","Calumenin")),symnum.args = signif_labels,aes(label = paste0("p = ", after_stat(p.format))),
                     label.x.npc="middle",label.y.npc="top") +facet_grid()
Ag_plot

file_name="/ix/djishnu/Trirupa/Schisto_Proj/ER_run/SEAplus_Eggneg/Plots/unsup_class2_ER/version2_30May23/SEA_SM25.Agmodel.LASSO_rep20_EndoA_B.pdf"
ggsave(filename=file_name,plot=Ag_plot,device="pdf",dpi=300,width=4, height=4,limitsize = FALSE)



