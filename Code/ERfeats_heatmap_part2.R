### This code plots heatmaps 

#rm(list = ls())
#cat("\014")

library(ggplot2) # ggplot() for plotting
library(dplyr)
library(plyr)# data reformatting
library(tidyr) # data reformatting
library(stringr) # string manipulation
library(RColorBrewer)
library(sjPlot)
library(pheatmap)
library(genefilter)
###creating features X patId matrix
patID_Xdf=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/forHeatmap/data/c2A_2B_featsorted_30May23_V2_patid_X.csv", row.names=1)
patID_X=rev(patID_Xdf)
#patID_X=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/forHeatmap/data/c2A_2B_featsorted_30May23_V2_patid_X.csv", row.names=1)
#colnames(patID_X) <- gsub('.', '-', colnames(patID_X))
#all_feats=c(feat_Z1,feat_Z2,feat_Z3,feat_Z4,feat_Z5)
all_feats=colnames(patID_X)

#patID_X_feat <- patID_X[, which((names(patID_X) %in% all_feats)==TRUE)]

X_feat_patID=as.data.frame(t(patID_X))
#X_feat_patID=as.data.frame(t(patID_X_feat)) ##this is the tmpMat
tpmMat=X_feat_patID
##creating patID X group labels matrix. This is colAnnot
## sort the features in the feature matrix based on heatmap_featMFIsort.py script and then use the input ##
patID_grp=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/forHeatmap/data/c2A_2B_featsorted_30May23_V2_patid_Y.csv", row.names=1)
feat_list=data.frame(target_id=all_feats,feat_name=all_feats)
#gene_list=feat_list

results <- list()
# Z-transform
tpmMat <- as.matrix(tpmMat)
row_means <- rowMeans(tpmMat)
row_stds <- rowSds(tpmMat)
zMatrix <- (tpmMat - row_means) / row_stds

colnames(zMatrix)=paste0("id_",colnames(zMatrix))
rownames(patID_grp)=paste0("id_",rownames(patID_grp))
# subset z-matrix for selected feats
subZMat <- data.frame(zMatrix[rownames(zMatrix) %in% feat_list$target_id, ])
## browser()
subZMat_ordered <- cbind(subZMat[feat_list$target_id, ], "feat_name"= feat_list$feat_name)

rownames(subZMat_ordered) <- subZMat_ordered$feat_name
subZMat2 <- subZMat_ordered[,-ncol(subZMat_ordered)] ##
results[["zMat"]] <- subZMat2

# Define color scale
margins <- c(min(subZMat2), max(subZMat2)) 
breaksList=NA
if(sum(abs(margins) > 1.5) > 0){
  breaksList <- seq(-1.5, 1.5, by = 0.03)
} 

feat_list=data.frame(target_id=rownames(subZMat2),feat_name=rownames(subZMat2))
#gene_list=feat_list

patID_grp[patID_grp$Y == 1, "Y"] <- "EndoA"
patID_grp[patID_grp$Y == 0, "Y"] <- "EndoB"
patID_grp$Y=as.factor(patID_grp$Y)
colAnnot=patID_grp

##function for saving heatmap
save.heatmap <- function(x, filename,width=6,height=2.5) {
  pdf(filename,width=width,height=height) 
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


# heatmap
savePath="/ix/djishnu/Trirupa/Schisto_Proj/forHeatmap/plots/c2A_c2B.ScAg.MFI50._spec0.1_heatmapSorted_v1.pdf"
annot_colors=list(Y=c("EndoA"="#3B9AB2","EndoB"="#FFBF00"))

hMap <- pheatmap(as.matrix(subZMat2), color=colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100), 
                 breaks=breaksList,show_rownames=TRUE, show_colnames=FALSE,
                 cluster_rows = FALSE, rownames=rownames(subZMat2),cluster_cols = FALSE,
                 annotation_col = colAnnot, fontsize_row = 4,fontsize_col = 4,annotation_colors=annot_colors,annotation_legend=TRUE,
                 legend=TRUE, border_color = "grey", fontsize = 6)


save.heatmap(x = hMap, filename = savePath)
results[["hMap"]] <- hMap




