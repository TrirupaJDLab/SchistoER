## This code plots the MFI values for the endotypes and SEA+Egg+ ###

#rm(list = ls())
#cat("\014")

library(dplyr)
library(pROC)
library(ROCR)
library(ggplot2)
library(reshape2)
library(ggpubr)

xydf=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/data/c3_c2lab0lab1_AbScAg.MFI50_30May23_XY_patid.V2.csv")

##create lists 
Z1=c("IgG.Gst.SM25", "IgG1.Gst.SM25", "IgG2.SEA", "IgG3.MEG", "IgG3.CD63", "FcR1.Gst.SM25", "FcR2A.Gst.SM25", "FcRIIB.Gst.SM25", "FcRIIIA.Gst.SM25", "FcR3B.Gst.SM25", "SNA.SEA")
Z2=c("IgG.MEG", "IgG.Calumenin", "IgG2.Gst.SM25", "IgG3.CD63", "IgG4.SEA", "IgG4.CD63", "FcREpsilonI.SEA", "FcR1.SEA", "FcR2A.SEA", "FcR2A.Calumenin", "FcR3B.Calumenin")
Z3=c("IgG.Gst.SM25", "IgG1.SEA", "IgG1.Gst.SM25", "IgA2.Gst.SM25", "FcR1.SEA", "FcR1.Gst.SM25", "FcR2A.Gst.SM25", "FcRIIB.Gst.SM25", "FcRIIIA.Gst.SM25", "FcR3B.Gst.SM25", "SNA.SEA", "SNA.Gst.SM25", "SNA.CD63")
Z4=c("IgG.Gst.SM25", "IgG1.Gst.SM25", "IgG4.SEA", "IgG4.CD63", "FcREpsilonI.SEA", "FcRIIB.Gst.SM25")
Z5=c("IgG.Calumenin", "IgG1.Calumenin", "FcR1.CD63", "FcR1.Calumenin", "FcR2A.Calumenin", "FcRIIB.Calumenin", "FcRIIIA.Calumenin", "FcR3B.Calumenin")

allZfeats=c(Z1,Z2,Z3,Z4,Z5)
ERfeats=unique(allZfeats)


subset_XY=xydf %>% select(all_of(ERfeats))
subset_XY=log(subset_XY,2)
subset_XY$Y=xydf$Y
##renaming Ys
subset_XY$Y[subset_XY$Y == 2] <- "SEA+ Egg+"
subset_XY$Y[subset_XY$Y == 1] <- "Endo 2A"
subset_XY$Y[subset_XY$Y == 0] <- "Endo 2B"

signif_labels=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

#subset_melt=melt(subset_XY)
my_comparisons <- list( c("Endo 2A", "Endo 2B"), c("Endo 2B", "SEA+ Egg+"), c("Endo 2A", "SEA+ Egg+"))
annot_colors=c("SEA+ Egg+"="#E15759","Endo 2A"="#3B9AB2", "Endo 2B"="#FFBF00")
Xi=data.frame()
for (i in ERfeats){
  Xi=subset_XY[,c("Y",i)]
  melt_Xi=melt(Xi)
  #colnames(melt_Xi)[3] <- "log2 MFI)"
  ##plotting the melted df for each i (feature)
  melt_plot=ggplot2::ggplot(data=melt_Xi,aes(x=Y,y=value))+
    ggplot2::geom_boxplot(data=melt_Xi, width=0.40, aes(fill=as.factor(Y))) + 
    scale_fill_manual(values=annot_colors) + 
    scale_x_discrete(labels=c('Endo 2A', 'Endo 2B', 'SEA+ Egg+')) + 
    ggtitle(as.character(i)) +
    theme(axis.text = element_text(size = 13),plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.spacing = unit(0.5, "lines"),
          panel.border = element_rect(fill = NA, color = "black",linetype="dashed"),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 15)) +
    ylab("log2(MFI)") + xlab("Status") +
    stat_compare_means(comparisons = my_comparisons, symnum.args = signif_labels,
                       aes(label = paste0("p =", after_stat(p.format))),label.x.npc="middle",label.y.npc="top") 
   ##### saving plot
   file_name=paste0("/ix/djishnu/Trirupa/Schisto_Proj/ER_run/SEAplus_Eggneg/Plots/unsup_class2_ER/version2_30May23/R_boxplots_3way/",i,"_boxplot_c2AB_c3_30May23.pdf")
   ggsave(filename=file_name,plot=melt_plot,device="pdf",dpi=300,width=5, height=6)
}

melt_plot


