### This code creats plots to compare the AUCs from randomly sampled LFs (equal to the number of actual LFs selected by ER) at different spec parameters to select the appropriate spec)

#rm(list = ls())
#cat("\014")

library(SLIDE)
library(pROC)

################################################################################
Z <- read.table("/ix/djishnu/Trirupa/Schisto_Proj/ER_run/Egg_NoEgg_revisedLabels/MFI50_filt/_30may23_f0.5ERpipe3_del0.1_lam1.output/zmatrix_EvsNE.MFI50.d0.1_lam1.auc.csv", sep = ",", header = T, row.names = 1)
y <- read.table("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/EvsNE_Ab_ScAg_30May23_Y.csv", row.names = 1, sep = ",", header = T)
colnames(y) <- "y"
zz <- row.names(Z)
y <- y[zz,,drop=F]
source("/ix/djishnu/Javad/Hierarchical_ER_v2/R/pairwise_interactions.R")


# Real ##########################################################################

# put your marginal and interaction variables here!
sigK <- c(4,5,9,13,16)
sigIn=c()
IntData <- pairwise_interactions(sigK,Z)
Dataint <- IntData$interaction[, sigIn]
Data_real <- data.frame(y = y, Z[, sigK], Dataint)
lmod  <- lm(y~.,data=Data_real)
yhat <- predict(lmod,Data_real[,-1],type = 'response')
aucreal <- auc(response=as.matrix(y), predictor=as.matrix(yhat))
aucreal_2=0.9361

# All random ####################################################################
Fullreport <- NULL

for (i in 1:1000) {
  sigKRandom <- sample(ncol(Z), size = length(sigK)) ## Random marginal
  Data_fullRandom <- data.frame(y = y, Z[, sigKRandom])
  SumInt <- summary(lm(y ~ ., data = Data_fullRandom))
  lmod  <- lm(y~.,data=Data_fullRandom)
  yhat <- predict(lmod,Data_fullRandom[,-1],type = 'response')
  aucrandom <- auc(response=as.matrix(y), predictor=as.matrix(yhat))
  Fullreport <- rbind(Fullreport, aucrandom)
}


################################################################################
## Report
library(ggplot2)
library(reshape2)
rawdf <- data.frame(FullRandom = Fullreport)
df <- melt(rawdf)
colnames(df) <- c("group", "value")

cols <- c("#0000FF", "#00FF00")

# Basic density plot in ggplot2

P2 <- ggplot(df, aes(x           = value, fill = group)) +
  geom_density(alpha            = 0.7) +
  scale_fill_manual(values      = cols) +
  theme_light() +
  geom_vline(xintercept         =aucreal , 
             linetype = "longdash",
             colour = "red",size=2) + xlab("auc")

P2 <- P2 + annotate("text", x     = aucreal + 0.01, 
                    y = 55, 
                    label = " ", 
                    angle = "90") +
  xlim(0.25, 1.0)
P2 <- P2 + theme(panel.border = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 axis.text = element_text(size = 20),
                 axis.title.x =element_text(size = 20),
                 axis.title.y = element_text(size = 20))




P2

saveRDS(df, file                = "/ix/djishnu/Trirupa/Schisto_Proj/ER_run/Egg_NoEgg_revisedLabels/MFI50_filt/_30may23_f0.5ERpipe3_del0.1_lam1.output/ERplots/ggplotobject_randompar_data2.rds")
ggsave(P2, filename              = "/ix/djishnu/Trirupa/Schisto_Proj/ER_run/Egg_NoEgg_revisedLabels/MFI50_filt/_30may23_f0.5ERpipe3_del0.1_lam1.output/ERplots/randomvsreal_specPlot.pdf")

