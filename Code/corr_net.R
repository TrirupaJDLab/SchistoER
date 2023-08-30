### this code makes correlation networks for each LF with input as the nodes of the LF and the original feature matrix ###

####
#install.packages("qgraph")
# load packages
library(ggplot2) # ggplot() for plotting
library(dplyr) # data reformatting
library(tidyr) # data reformatting
library(stringr)

xdf=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/EvsNE_AbScAg_30May23_MFI50_X.csv",row.names=1)
#scaledData <- scale(xdf[,-1])
lis=xdf[,c("IgG.Gst.SM25", "IgG1.Gst.SM25" , "IgG2.SEA", "IgG3.MEG")]
library(qgraph)
# create example correlation matrix
corr_matrix <- cor(lis)
# create correlation network using qgraph
corr_plot=qgraph(corr_matrix, filename="corrnet_Z4_v3_thres0.5",layout = "spring", threshold=0.4,repulsion=0.4,labels = colnames(corr_matrix),label.font=2,label.scale.equal=TRUE,label.prop=0.95,shape="ellipse", posCol="#4E79E8", negCol="#59A14F",filetype='pdf',height=5,width=10)

corr_plot
