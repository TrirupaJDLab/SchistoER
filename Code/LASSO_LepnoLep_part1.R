####
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

dataset=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/data/3way_AbScAg_30May23_MFI50.csv") ##has both X (outcome) and Y (labels)
data=dataset[,-1]
Y_data <- data[,1]
colnames(data)[1] <- "Y" 
X_data=data[,-1]

################Permutation Test Outside##############

################ Set K for K-fold CV  #######################
k <- 10
aucRF <- c()
aucRF_permut <- c()
aucSVM <- c()
aucSVM_permut <- c()
##lists for saving aucs from each rep
aucRF_rep = c()
aucSVM_rep = c()
aucRF_permut_rep = c()
aucSVM_permut_rep = c()
##lists for saving AUCs and features
trueY <- c()
permutY <- c()

append.RF <- c()
append.SVM <- c()

append.RF_permut <- c()
append.SVM_permut <- c()

#selectedVars <- c()


dF <- c(1:nrow(data))
permutedY <- sample(Y_data) ##sampling Y from the main dataset

num_repCV=20
#replicateCV=1
for (replicateCV in 1:num_repCV){
  #print("in Replication")
  ################### Scramble The Response Vector ###############
  permutedY <- sample(Y_data) 
  
  ################## Get Folds for the Data ######################
  folds <- createFolds(y = data$Y, k=k, list = FALSE, returnTrain = FALSE)
  # at each fold, we make sure to start with the true y label
  ori_myData <- cbind(data, folds)
  

  
  #for each fold of the data
  #NoF=1
  for (NoF in 1:k){
    fold = which(folds == NoF)
    #print(dim(ori_myData))
    myData <- ori_myData
    
    
    #1 is regular, 2 is permute 
    for (m in 1:2){
      # we set y as regular
      if (m == 1){
        #print('Regular Lasso Running...')
        train <- myData[-fold, ]
        test <- myData[fold, ]
        
        # delete the fold column
        train <- train[, -ncol(train)]
        test <- test[, -ncol(test)]
        
        X = as.matrix(train[, -1]) # X stays the same for permut and non-permut
        y_train = train$Y
      }
      # the second loop, we use scrambled y
      if (m == 2){
        #print('Permutation Test Lasso Running...')
        #print(dim(myData))
        myData$Y <- permutedY
        train <- myData[-fold, ]
        test <- myData[fold, ]
        
        # delete the fold column
        train <- train[, -ncol(train)]
        test <- test[, -ncol(test)]
        
        X = as.matrix(train[, -1]) # X stays the same for permute and non-permute
        y_train = train$Y
      }
      
      
      #feature selection
      variables <- c()
      #print("selecting features")
      ########################### Loop for Feature Selection ###############
      #num_FS=1
      for (num_FS in 1:20){
        #print(num_FS)
        ######################### Run CV.Glmnet #############################
        glmnet1 <- cv.glmnet(X, y=train$Y, alpha=1, nfolds = k)
        
        # tuning lambda
        lambda <- glmnet1$lambda.min
        #lambda <- lambda*tun_factor
        lambda <- lambda*0.5
        
        ###################### Using Selected Lambda to Run Lasso ############
        glmnet2 <- glmnet(X, y=train$Y, alpha=1, lambda = lambda)
        c <- coef(glmnet2)
        
        ##################### Saving Selected Features From Lasso ############
        inds<-which(c!=0)
        tmp_variables <- tail(row.names(c)[inds], -1)
        len=length(tmp_variables)
        
        #selectedVars <- append(selectedVars, len)
        #remove intercept from the list
        # tmp_variables <- tail(row.names(c)[inds], -1)
        #print(tmp_variables)
        ## if ....
        variables<-c(variables,tmp_variables)

      }
      ################# Selecting Features Based On Frequency ############
      freq = sort(table(variables),decreasing=TRUE)/num_FS
      #print(freq)
      
      stable_Var = names(which(freq>0.6))  #### the 0.6 can be changed ######
      cat(replicateCV,"NoF:",NoF,length(stable_Var))
      print(paste0("m=",m))
      print(stable_Var)
      
      #selectedVars <- append(selectedVars, length(stable_Var))
      ###selectedVars <- append(selectedVars, stable_Var)
      ################## Use Stable Features on RF & SVM ##################### 
      tempTr <- train[, (names(train) %in% stable_Var)]
      
      
      tempTrain <- cbind.data.frame(train$Y, tempTr)
      
      
      tempTr <- test[, (names(test) %in% stable_Var)]
      
      
      tempTest <- cbind.data.frame(test$Y, tempTr)
      
      
      colnames(tempTrain)[1] <- "Y"
      colnames(tempTest)[1] <- "Y"
      
      ######################## Run SVM and RF #########################
      # append all test Y for AUC calculation   
      #Y <- append(Y, tempTest$Y)
      #For true Y and permut Y, we put them in separate lists but name both Y at the end for readability
      if (m == 1){
        trueY <- append(trueY, tempTest$Y)
        Y_train <- train$Y ##this Y will be used to train SVM and RF
      }
      if (m == 2){
        permutY <- append(permutY, tempTest$Y)
        Y_train <- train$Y ##this Y will be used to train SVM and RF
      }
      
      #svm
      colnames(tempTrain)[1]<- "Y_train"
      svmfit = svm(Y_train~ ., data = tempTrain , kernel="linear", cost=19, scale=FALSE, tolerance=0.01)
      yhat.SVM = predict(svmfit, newdata = tempTest[, -1])
      
      #rf
      RFfit <- randomForest(Y_train~ ., data = tempTrain, importance=TRUE, ntree = 15)
      yhat.RF = predict(RFfit, newdata = tempTest[, -1])
      
      #################### Store yhat in different array depend on the mode #############
      if (m == 1){
        append.SVM <- append(append.SVM, yhat.SVM)
        append.RF <- append(append.RF, yhat.RF)
      }
      if (m == 2){
        append.SVM_permut <- append(append.SVM_permut, yhat.SVM)
        append.RF_permut <- append(append.RF_permut, yhat.RF)
      }
    }
  }
  # get auc score for each replication
  #print("printing AUCs")
  aucRF_rep[replicateCV] <- auc(trueY, append.RF)
  #print(aucRF[replicateCV])
  aucSVM_rep[replicateCV] <- auc(trueY, append.SVM)
  #print(aucSVM[replicateCV])
  
  aucRF_permut_rep[replicateCV] <- auc(permutY, append.RF_permut)
  #print(aucRF_permut[replicateCV])
  aucSVM_permut_rep[replicateCV] <- auc(permutY, append.SVM_permut)
}
# get auc score for each replication
aucRF <- auc(trueY, append.RF)
#print(aucRF[replicateCV])
aucSVM<- auc(trueY, append.SVM)
#print(aucSVM[replicateCV])

aucRF_permut <- auc(permutY, append.RF_permut)
#print(aucRF_permut[replicateCV])
aucSVM_permut <- auc(permutY, append.SVM_permut)

write.csv(aucRF_rep, "/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/aucRF_allpat_AbScAg_MFI50_l0.5.csv")
write.csv(aucSVM_rep, "/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/aucSVM_allpat_AbScAg_MFI50_l0.5.csv")
write.csv(aucRF_permut_rep, "/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/aucRFpermute_allpat_AbScAg_MFI50_l0.5.csv")
write.csv(aucSVM_permut_rep, "/ix/djishnu/Trirupa/Schisto_Proj/LASSO/AUCs/EvsNE_lepNolep/aucSVMpermute_allpat_AbScAg_MFI50_l0.5.csv")