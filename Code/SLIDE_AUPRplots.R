## loading libraries 
library(PRROC)
library(ggplot2)
library(reshape2)
library(ggpubr)


analysis_list=c("EggvsNoEgg","class2_active","Endotypes")
for (analysis in analysis_list){
  yaml_path=paste0("./Code/pipeline3_",analysis,".yaml")

  er_input <- yaml::yaml.load_file(yaml_path)

  pathLists <- list.files(er_input$out_path,recursive = T,pattern = "results")
  perfList <- lapply(paste0(er_input$out_path,pathLists), readRDS)

  AUPR_scoresDF = data.frame(nrep = numeric(), actual_PR = numeric(), perm_PR = numeric())

  for (i in 1:length(perfList)){
    fold_df=as.data.frame(perfList[[i]][["each_fold"]])
    fold.ER_DF=fold_df[fold_df$method == c("plainER","plainER_y"),]
  
    # Extract pred_vals where method is "plain_ER"
    actual_pred_vals <- fold.ER_DF$pred_vals[fold.ER_DF$method == "plainER"]
    actual_true_vals <- fold.ER_DF$true_vals[fold.ER_DF$method == "plainER"]
  
    perm_pred_vals <- fold.ER_DF$pred_vals[fold.ER_DF$method == "plainER_y"]
    perm_true_vals <- fold.ER_DF$true_vals[fold.ER_DF$method == "plainER_y"]
  
    # Calculate Precision-Recall curve
    actual_pr_curve <- pr.curve(scores.class0 = actual_pred_vals, weights.class0 = actual_true_vals, curve = TRUE)
    actual_PR = actual_pr_curve$auc.integral
    
    perm_pr_curve <- pr.curve(scores.class0 = perm_pred_vals, weights.class0 = perm_true_vals, curve = TRUE)
    perm_PR = perm_pr_curve$auc.integral
    
    ## updating values in the df
    append_vals = c(i,actual_PR,perm_PR)
    AUPR_scoresDF = rbind(AUPR_scoresDF,append_vals)
  }

  colnames(AUPR_scoresDF)=c("nrep","actual_PR","perm_PR")
  melted_DF= melt(AUPR_scoresDF, measure.vars = c("actual_PR", "perm_PR"))
  colnames(melted_DF)=c("nrep","Method","AUPR")
  
  signif_labels=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  lambda_boxplot <- ggplot2::ggplot(data = melted_DF,
                                    ggplot2::aes(x = Method, y = AUPR, fill = Method)) +
    ggplot2::geom_boxplot(width = 0.20) +
    ggplot2::labs(fill = "Method") + scale_x_discrete(labels = c('Actual', 'Permuted')) +
    ggplot2::scale_alpha(guide = 'none') +
    ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    stat_compare_means(symnum.args = signif_labels,
                       aes(label = paste0("p = ", after_stat(p.format))),
                       label.x.npc = "middle", label.y.npc = "top") +
    ggplot2::ggtitle(analysis) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  lambda_boxplot
  
  ggsave(paste0("ER_AUPR_",analysis,".pdf"), plot = lambda_boxplot, width = 6, height = 4, units = "in")
}
