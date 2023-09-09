# SchistoER
#The Code folder contains all the scripts and yaml files used to perform the analysis.

##The Input folder contains the following folders:

1. Validation-  contains the inputs for performing ER crossprediction, including the training dataset (Brazilian discovery cohort) and the test dataset (Kenyan validation cohort).
2. X_data- has processed .csv files with the feature the inputs for ER (EvsNE : Egg+ vs Egg- ; class2: SEA+Egg-)
3. Y_data- Has the outcome labels for the corresponding X datasets, as inputs for the ER pipeline. (class2: SEA+Egg-)
4. unsupervised_XYdata- has inputs with labels generated from unsupervised clustering (k-means)

##The Code folder contains all scripts used for the analyses:

1. Scripts for the ER pipeline:
i) the files ending in ".yaml", are to be tuned (mainly for delta, lambda values) for each analysis.
ii) pipe3_Schisto.R : the R script for running ER. Has .yaml files as the inputs and the outputs for each analysis are in the output folder.

3wayunivariate_boxplot.R : script to plot the univariate differences in MFI for EndoA, EndoB and SEA+Egg+ (active) samples. It takes as input the MFI values for each feature along with the outcome labels, with the specified significant latent factors (Z) obtained from ER (the features selected in each).

ER_crossprediction.Rmd : has the R markdown file for ER crossprediction between the discovery (Brazilian) cohort and the validation (Kenyan) cohort. 
it takes the following inputs:
a) X_train_AllERselected_featEvsNE_30May23.csv : file containing MFI values of features selected by the significant latent factors in the discovery cohort.
b) Y_train_AllERselected_featEvsNE_30May23.csv: file with outcome labels (Y) in the discovery cohort.
c) XY_ScAg.VC_EvsNE.data_imputed_19jun23.csv: file containing the MFI values of the features in the validation cohort within the significant latent factors (obtained from discovery cohort) and the Y labels.
d) ER output: the rds file generated by running pipe3_Schisto.R
e) significant Latent factors selected by running ER on the discovery cohort.

ERfeats_directionality.R : This script computes AUCs of features picked up by the significant latent factors, for both discovery and validation cohorts. Inputs are, X and Y .csv files for the respective analysis in the discovery cohort, along with the dataset having X and Y from the validation cohort, and the significant latent factors from the discovery cohort.

ER_heatmap_part1.py : this script orders the features selected by the significant latent factors based on the mean of log2(fold change) between the 2 groups (Y) being compared. 

ER_heatmap_part2.R : this script takes the output of ER_heatmap_part1.py and plots a heatmap. The inputs of this script are the following:
a) c2A_2B_featsorted_30May23_V2_patid_X.csv: This file has patient ids and the corresponding feature MFIs. 
b) c2A_2B_featsorted_30May23_V2_patid_Y.csv: This file has only the patient ids and the corresponding outcome (labels).

ERspec_controlPlot.R : 
This code creats plots to compare the AUCs from randomly sampled latent factors (equal to the number of actual latent factors selected by ER, at different spec parameters to select the appropriate spec). This code takes as inputs the matrix of latent factors obtained from ER (in script XXXX) along with the Y (outcome labels) for the corresponding analysis.

LASSO_Ag_models_withCV.R : this code performs LASSO regression using antigens highlighted by the significant latent factors from ER. It builds and compares antigen-specific models in predicting Schisto outcome (Y). The inputs for this script includes X (feature dataset) and Y (outcome/labels) for the corresponding analysis.

Unsup_clustering.py : This code does unsupervised clustering by computing a silhouette index to identify the optimal k , and then generates labels using k-means clustering. The input required is the .csv file containing class2 (SEA+Egg-) samples i.e., the feature set (X) and the outcome Y, (labels). 

ZLevel_boxplot.R : This script plots a boxplot illustrating the distributions of the outcome labels by each significant latent factor. The following are inputs required to run the script:
1) x: .csv file with the MFI for features corresponding to the comparison of interest.
2) y: .csv file having the outcome (labels) for the comparisons (in binary; Egg+ denoted by 1, Egg- denoted by 0)
3) z: .csv file containing the matrix of all latent factors obtained from the ER pipeline.

corr_net.R : this code makes correlation networks for features within a latent factor. The inputs include:
1) the dataframe of feature MFI values
2) list of features selected by the significant latent factor. 

