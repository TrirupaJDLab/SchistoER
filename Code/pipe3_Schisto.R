### This code runs the ER pipeline with input as the yaml files for corresponding analysis ###
library(EssReg)
library(doParallel)
library(dplyr)
library(pROC)
library(ROCR)

registerDoParallel(detectCores())


###yaml path for check_runs
yaml_path = "pipeline3_EggvsNoEgg.yaml"
pipe3_run <- pipelineER3(yaml_path)
pipe3_run
