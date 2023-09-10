### this code makes correlation networks for each LF with input as the nodes of the LF and the original feature matrix ###

####
#install.packages("qgraph")
# load packages
library(ggplot2) 
library(dplyr) 
library(tidyr) 
library(stringr)
library(qgraph)

xdf=read.csv("EvsNE_X.csv",row.names=1)
lis=xdf[,c("IgG.Gst.SM25", "IgG1.Gst.SM25" , "IgG2.SEA", "IgG3.MEG")]  ### insert features in a particular significant LF

# create correlation matrix
corr_matrix <- cor(lis)
# create correlation network using qgraph
corr_plot=qgraph(corr_matrix, filename="corrnet_Z4",layout = "spring", threshold=0.4,repulsion=0.4,labels = colnames(corr_matrix),label.font=2,label.scale.equal=TRUE,label.prop=0.95,shape="ellipse", posCol="#4E79E8", negCol="#59A14F",filetype='pdf',height=5,width=10)

corr_plot
