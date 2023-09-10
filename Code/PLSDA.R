### this code makes PLS-DA plots for features selected in significant latent factors from ER ####

library(pls)
library(ggplot2)
library(pROC)
library(cvAUC)

xdf=read.csv("class2_active_X.csv", row.names=1)
x_std <- scale(xdf, T, T)
ydf=read.csv("class2_active_X.csv", row.names=1)
Dataset = cbind(ydf,x_std)

## significant Zs from ER pipeline
Z1=c("IgG1.SEA", "IgG1.Gst.SM25", "IgG3.SEA", "IgG3.Calumenin", "IgA.Gst.SM25", "IgA1.SEA", "IgA1.MEG", "IgA2.MEG", "IgA2.Calumenin", "FcR2A.Calumenin","FcR3B.Calumenin", "SNA.MEG", "SNA.CD63")
Z2=c("IgG2.Gst.SM25", "IgG2.Calumenin", "IgG3.CD63", "IgG4.SEA", "IgG4.CD63", "FcREpsilonI.SEA", "FcR2A.Calumenin")
Z3=c("IgG.Calumenin", "IgA.Calumenin", "IgA1.SEA", "IgA1.Gst.SM25", "IgA1.MEG", "IgA1.CD63", "IgA1.Calumenin", "IgA2.Calumenin")
Z4=c("IgG.Calumenin", "IgG1.Calumenin", "FcR1.CD63", "FcR1.Calumenin", "FcR3B.Calumenin")
allZfeats=c(Z1,Z2,Z3,Z4)
features=unique(allZfeats)

Y <- Dataset$Y
Outcome <- Dataset$Y

loc1 <- which(Outcome == 1)
loc0 <- which(Outcome == 0)

Outcome[loc1] <- "SEA+ Egg+"
Outcome[loc0] <- "SEA+ Egg-"

data <- Dataset[, features]

# Partial Least Squares

data <- cbind.data.frame(Y, data)
pls.fit = plsr(Y~., data = data, scale=TRUE)
summary(pls.fit)
pls.fit$loadings

AUC(pls.fit$scores[,1], Dataset$Y)
multiclass.roc(pls.fit$scores[,1], Dataset$Y)


PC1 <- pls.fit$scores[,1]
PC2 <- pls.fit$scores[,2]
dF <- cbind.data.frame(Outcome, PC1, PC2)
dF$Outcome <- as.factor(dF$Outcome)
p <- ggplot(dF, aes(x=PC1, y=PC2, colour = Outcome)) + geom_point(size = 2.5)
p +theme_classic()+ scale_color_manual(values=c("#EDC948","#E15759"))+ xlab("LV1") + ylab("LV2")

ggsave("PLS_class2_active.svg")
