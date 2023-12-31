# SchistoAbomics

## Input 
NOTE: All input datasets provided are pre-filtered for features based on median MFI>50. 

1. Validation-  contains the inputs for performing ER crossprediction, including the training dataset (Brazilian discovery cohort) and the test dataset (Kenyan validation cohort).
2. X_data- has processed .csv files with the feature the inputs for ER (EvsNE : Egg+ vs Egg- ; class2: SEA+Egg-)
3. Y_data- Has the outcome labels for the corresponding X datasets, as inputs for the ER pipeline. (class2: SEA+Egg-)
4. unsupervised_XYdata- has inputs with labels generated from unsupervised clustering (k-means)
5. Leprosy_confounding- has .csv files with MFI values and Leprosy and Schisto status. Allpat.csv includes all patients, PatNoLep.csv includes patients without leprosy. aucRFpermute*.csv files contain the permuted AUCs obtained from LASSO (with Random forest) for all patients and with patients without leprosy. aucRF_*.csv files contain AUC values for actual labels for the aforementioned 2 categories. 

## Code 
This folder contains all scripts and yaml files used to perform the analysis used for the analyses:

### 1. Scripts for the ER pipeline:
i) the files ending in ".yaml", are to be tuned (mainly for delta, lambda values) for each analysis.<br>
ii) pipe3_Schisto.R : the R script for running ER. Has .yaml files as the inputs and the outputs for each analysis are in the output folder.<br>
iii) ERmodel_performancePlot.R : makes box plots comparing AUCs on Actual and Permuted labels using the output of ER as input (pipeline_step5.rds).<br>
iv) LatentFactor_extraction.R : this code gets the signficant LFs and features within it from the output of ER pipeline. The inputs are X (feature MFIs) and Y (outcome/labels) of the comparison of interest; ER output (final_delta_0.1_lambda_1.rds).

### 2. ER_crossprediction.Rmd : 
This is an R markdown file for ER crossprediction between the discovery (Brazilian) cohort and the validation (Kenyan) cohort. It takes the following inputs:<br>
i) X_train.csv : file containing MFI values of features selected by the significant latent factors in the discovery cohort.<br>
ii) Y_train.csv : file with outcome labels (Y) in the discovery cohort.<br>
iii) XY_VC_imputed.csv: file containing the MFI values of the features in the validation cohort within the significant latent factors (obtained from discovery cohort) and the Y labels.<br>
iv) final_delta_0.1_lambda_1.rds: the rds file generated by running pipe3_Schisto.R on X_train.csv and Y_train.csv.<br>
v) significant Latent factors selected by running ER on the discovery cohort.

### 3. 3wayunivariate_boxplot.R : 
This script plots the univariate differences in MFI for EndoA, EndoB and SEA+Egg+ (active) samples. It takes as input the MFI values for each feature along with the outcome labels, with the specified significant latent factors (Z) obtained from ER (the features selected in each).

### 4. ERfeats_directionality.R : 
This script computes AUCs of features picked up by the significant latent factors, for both discovery and validation cohorts. Inputs are- X and Y .csv files for the respective analysis in the discovery cohort, along with the dataset having X and Y from the validation cohort, and the significant latent factors from the discovery cohort.

### 5. ER_heatmap_part1.py : 
This script orders the features selected by the significant latent factors based on the mean of log2(fold change) between the 2 groups (Y) being compared. 

### 6. ER_heatmap_part2.R : 
This script takes the output of ER_heatmap_part1.py and plots a heatmap. The inputs of this script are the following:<br>
i) endoA_endoB_featsorted_patid_X.csv: This file has patient ids and the corresponding feature MFIs. <br>
ii) endoA_endoB_featsorted_patid_Y.csv: This file has only the patient ids and the corresponding outcome (labels).

### 7. ERspec_controlPlot.R : 
This code creats plots to compare the AUCs from randomly sampled latent factors (equal to the number of actual latent factors selected by ER, at different spec parameters to select the appropriate spec). This code takes as inputs the matrix of latent factors obtained from ER (in script LatentFactor_extraction.R) along with the Y (outcome labels) for the corresponding analysis.

### 8. LASSO_Ag_models_withCV.R : 
This code performs LASSO regression using antigens highlighted by the significant latent factors from ER. It builds and compares antigen-specific models in predicting Schisto outcome (Y). The inputs for this script includes X (feature dataset) and Y (outcome/labels) for the corresponding analysis.

### 9. Unsup_clustering.py : 
This code does unsupervised clustering by computing a silhouette index to identify the optimal k , and then generates labels using k-means clustering. The input required is the .csv file containing class2 (SEA+Egg-) samples i.e., the feature set (X) and the outcome Y, (labels). 

### 10. ZLevel_boxplot.R : 
This script plots a boxplot illustrating the distributions of the outcome labels by each significant latent factor. The following are inputs required to run the script:<br>
i) x: .csv file with the MFI for features corresponding to the comparison of interest. <br>
ii) y: .csv file having the outcome (labels) for the comparisons (in binary; Egg+ denoted by 1, Egg- denoted by 0). <br>
iii) z: .csv file containing the matrix of all latent factors obtained from the ER pipeline.

### 11. corr_net.R : 
This code makes correlation networks for features within a latent factor. The inputs include:<br>
i) the dataframe of feature MFI values. <br> 
ii) list of features selected by the significant latent factor. 

### 12. LASSO_LepnoLep_part1.R : 
This script performs LASSO regression (with RF/SVM) and generates AUCs for actual and permuted Y (outcome) values. The input is a .csv file containing the actual Y labels for egg status and the leprosy status (in binary, Lep+ :1; Lep- : 0). The code generates permuted and actual AUCs in separate files. 

### 13. LASSO_LepnoLep_part2.R : 
This code makes box plots for the results obtained from LASSO_LepnoLep_part1.R. Takes as input .csv files having AUCs of actual and permuted, generated from part1.

### 14. PLSDA.R :
This code makes PLS-DA plots for features selected in significant latent factors from ER.

## Output
This folder contains outputs from the ER pipeline used for analyses described in the paper. 


