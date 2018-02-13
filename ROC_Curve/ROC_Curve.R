##### ROC
##### Qiwen Kang
##### 03/04/2017


# Description -------------------------------------------------------------
# Here we use the same data sets from RF distance part


#### Loading packages
library(ape)
library(phangorn)
library(kdetrees)
library(parallel)

setwd("/home/qiwen/r/working/rd10_dna")

roc_test <- function(threshold){

  nonOut_index <- sample(1000,100)
  out_index <- sample(1000,1)
  
  oriValue_kde <- matrix(c(rep(0,100),1,rep(0,101)),ncol=2)
  oriValue_cu <- matrix(c(rep(0,100),1,rep(0,101)),ncol=2)
  
  
  normTree_cu <- paste("data/depth",threshold,"/set1",sep = "")
  normTree_kde <- paste("depth",threshold,"/set1",sep = "")
  normTree <- list.files(normTree_kde, pattern = ".tre$")

  ## This part is for KDE tree
  tree <- rmtree(101,10)
  treeDir <- paste(normTree_kde, "/", normTree[nonOut_index], sep="")
  tree[1:100] <- lapply(treeDir,read.tree)
  
  tree[[101]] <- read.tree(paste(normTree_kde, "/", normTree[1000+out_index], sep=""))
  kdeResul <- kdetrees(tree, k = 1.5)
  oriValue_kde[kdeResul$i,2] <- 1
  
  # True positive
  tp_kde <- ifelse(oriValue_kde[101,2]==1, 1, 0)
  # False negative
  fn_kde <- ifelse(oriValue_kde[101,2]==0, 1, 0)
  # True positive rate
  tpr_kde <- tp_kde/(tp_kde + fn_kde)
  
  # True negative
  tn_kde <- 100 - sum(oriValue_kde[1:100,2])
  # False positive
  fp_kde <- sum(oriValue_kde[1:100,2])
  # False positive rate
  fpr_kde <- fp_kde/(fp_kde + tn_kde)
  
  
  ## This part is for CU 
  nonOut_index <- sample(1000,100)
  out_index <- sample(1000,1)
  nonOutlier <- read.table(paste(normTree_cu,"/All_1.0_tree1.txt",sep=""))
  outlier <- read.table(paste(normTree_cu,"/All_2.0_tree1.txt",sep=""))
  (nonOutlier_ratio <- nonOutlier[,2][nonOut_index])
  (outlier_ratio <- outlier[,2][out_index])
  
  sample_set <- append(nonOutlier_ratio,outlier_ratio)
  
  # Detect the outlier
  indexTF_cu <- sample_set > quantile(sample_set,0.95)
  oriValue_cu[indexTF_cu,2] <- 1
  
  # True positive
  tp_cu <- ifelse(oriValue_cu[101,2]==1, 1, 0)
  # False negative
  fn_cu <- ifelse(oriValue_cu[101,2]==0, 1, 0)
  # True positive rate
  tpr_cu <- tp_cu/(tp_cu + fn_cu)
  
  # True negative
  tn_cu <- 100 - sum(oriValue_cu[1:100,2])
  # False positive
  fp_cu <- sum(oriValue_cu[1:100,2])
  # False positive rate
  fpr_cu <- fp_cu/(fp_cu + tn_cu)
  
  # The order is 
  result <- c(tpr_kde, fpr_kde, tpr_cu, fpr_cu)
  return(result)
}
options(mc.cores=6)
library(ggplot2)
# roc_point <- matrix(0,nrow = 50,ncol = 44)
# interest_value <- seq(-2,3,by=0.5)

roc_point <- matrix(0,nrow = 50,ncol = 24)
interest_value <- c("14","16","18","20","40","60")
# "06","07","08","09","10","12",
for(i in 1:50) roc_point[i,] <- unlist(lapply(interest_value,roc_test))

valueMatrix <- matrix(colMeans(roc_point),ncol = 4,byrow = T)
roc_tprKDE <- valueMatrix[,1]
roc_fprKDE <- valueMatrix[,2]
roc_tprCU <- valueMatrix[,3]
roc_fprCU <- valueMatrix[,4]

roc_kde_plot <- cbind(as.data.frame(roc_tprKDE),as.data.frame(roc_fprKDE))
roc_cu_plot <- cbind(as.data.frame(roc_tprCU),as.data.frame(roc_fprCU))
roc_plot <- rbind(cbind(unname(roc_kde_plot),method = "KDE"),cbind(unname(roc_cu_plot),method = "CU"))
names(roc_plot) <- c("TPR","FPR","Method")

save(roc_plot,file="data/plots/qiwenrocsim200.Rda")
ggplot(data = roc_plot, aes(x = FPR, y = TPR,lty=Method))+ 
  geom_line() +
  ggtitle("ROC Curve") 
