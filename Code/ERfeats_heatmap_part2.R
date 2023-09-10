### This code plots heatmaps 

library(ggplot2) 
library(dplyr)
library(plyr)
library(tidyr) 
library(stringr) 
library(RColorBrewer)
library(sjPlot)
library(pheatmap)
library(genefilter)

###creating features X patId matrix
patID_Xdf=read.csv("endoA_endoB_featsorted_patid_X.csv", row.names=1)
patID_X=rev(patID_Xdf)
all_feats=colnames(patID_X)


X_feat_patID=as.data.frame(t(patID_X))
tpmMat=X_feat_patID
##creating patID X group labels matrix. This is colAnnot
## sort the features in the feature matrix based on heatmap_featMFIsort.py script and then use the input ##
patID_grp=read.csv("endoA_endoB_featsorted_patid_Y.csv", row.names=1)
feat_list=data.frame(target_id=all_feats,feat_name=all_feats)
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
savePath="EndoA_EndoB_spec0.1_heatmap.pdf"
annot_colors=list(Y=c("EndoA"="#3B9AB2","EndoB"="#FFBF00"))

hMap <- pheatmap(as.matrix(subZMat2), color=colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100), 
                 breaks=breaksList,show_rownames=TRUE, show_colnames=FALSE,
                 cluster_rows = FALSE, rownames=rownames(subZMat2),cluster_cols = FALSE,
                 annotation_col = colAnnot, fontsize_row = 4,fontsize_col = 4,annotation_colors=annot_colors,annotation_legend=TRUE,
                 legend=TRUE, border_color = "grey", fontsize = 6)


save.heatmap(x = hMap, filename = savePath)
results[["hMap"]] <- hMap




