### This code plots SLIDE_cv model performance from SLIDE output. ##

library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsignif)

rds_pathList=c("./Output/SLIDEcv_EvsNE/delta0.05lambda1_boxplot_data.rds",
            "./Output/SLIDEcv_Class2vsActive/delta0.1lambda1_boxplot_data.rds",
            "./Output/SLIDEcv_EndoAvsEndoB/delta0.1lambda1_boxplot_data.rds")

for (i in seq_along(rds_pathList)){
  print(i)
  rds_path=rds_pathList[i]
  bp_df=readRDS(rds_path)
  bp_ER=bp_df %>% filter(!(method %in% c("lasso","lasso_y")))  
  Actual=as.vector(bp_ER[which(bp_ER$method == "SLIDE"),2]) ## these are unused in the code. The code automatically does it in wilcox.test()
  Permuted=as.vector(bp_ER[which(bp_ER$method == "SLIDE_y"),2])
  signif_labels=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  lambda_boxplot=ggplot2::ggplot(data = bp_ER,
                                 ggplot2::aes(x = method,
                                              y = auc,fill=method)) +
    ggplot2::geom_boxplot(width=0.20) + scale_x_discrete(labels=c('Actual', 'Permuted')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    stat_compare_means(symnum.args = signif_labels,aes(label = paste0("p = ", after_stat(p.format))),label.x.npc="middle",label.y.npc="top")+
    ggplot2::scale_alpha(guide = 'none') 
  
  print(lambda_boxplot)
  
  if (i==1){
    savepath="./Output/SLIDEcv_EvsNE/SLIDEcv_boxplot.pdf"
    ggsave(savepath,plot=lambda_boxplot,width=4,height=4, units="in")}
  if (i==2){
    savepath="./Output/SLIDEcv_Class2vsActive/SLIDEcv_boxplot.pdf"
    ggsave(savepath,plot=lambda_boxplot,width=4,height=4, units="in")}
  if (i==3){
    savepath="./Output/SLIDEcv_EndoAvsEndoB/SLIDEcv_boxplot.pdf"
    ggsave(savepath,plot=lambda_boxplot,width=4,height=4, units="in")}
 
}