##### GG plot for median 
##### Qiwen Kang
##### 01/18/2018

#### Set the workstation
setwd("/home/qiwen/r/180108SmallRFDistance/")

#### Read data sets.
cValue<-c("06","08","10","12","14","16","18","20","40","60")
sNum<-c(0:9)

fNum<-function(cValue,sNum){
  workDir<-paste("data/rd8/depth",cValue,"/set",sNum,"/",sep = "")
  
  fileList<-list.files(workDir,"^All")#
  geneOneTreeOne <- as.data.frame(read.table(paste(workDir,fileList[1],sep=""))[,2])
  geneOneTreeTwo <- as.data.frame(read.table(paste(workDir,fileList[2],sep=""))[,2])
  geneTwoTreeOne <- as.data.frame(read.table(paste(workDir,fileList[3],sep=""))[,2])
  geneTwoTreeTwo <- as.data.frame(read.table(paste(workDir,fileList[4],sep=""))[,2])
  
  
  ####Under Tree One
  
  fOneOne<-fivenum(geneOneTreeOne[,1])
  fTwoOne<-fivenum(geneTwoTreeOne[,1])
  fOneTwo<-fivenum(geneOneTreeTwo[,1])
  fTwoTwo<-fivenum(geneTwoTreeTwo[,1])
  treeOne<-ks.test(geneOneTreeOne[,1],geneTwoTreeOne[,1])
  treeTwo<-ks.test(geneOneTreeTwo[,1],geneTwoTreeTwo[,1])
  
  return(list("GeneOneTreeOne" = fOneOne, "GeneTwoTreeOne" = fTwoOne,
              "GeneOneTreeTwo" = fOneTwo, "GeneTwoTreeTwo" = fTwoTwo,
              "KS Test Tree One" = treeOne,
              "KS Test Tree Two" = treeTwo))
}
dataOneTreeOne<-list()
dataTwoTreeOne<-list()
dataTwoTreeTwo<-list()
dataOneTreeTwo<-list()
for(i in 1:10){
  
  res_06 <- fNum(cValue[1],sNum[i])
  res_08 <- fNum(cValue[2],sNum[i])
  res_10 <- fNum(cValue[3],sNum[i])
  res_12 <- fNum(cValue[4],sNum[i])
  res_14 <- fNum(cValue[5],sNum[i])
  res_16 <- fNum(cValue[6],sNum[i])
  res_18 <- fNum(cValue[7],sNum[i])
  res_20 <- fNum(cValue[8],sNum[i])
  res_40 <- fNum(cValue[9],sNum[i])
  res_60 <- fNum(cValue[10],sNum[i])
  
  #### data frame
  ### We set every set's third value since, in out case, we just care about the median
  ### 1:minumum 2:lower_hinge 3:median 4:upper-hinge 5:maximum
  dataOneTreeOne[[i]] <- rbind(res_06[[1]][3],res_08[[1]][3],
                               res_10[[1]][3],res_12[[1]][3],
                               res_14[[1]][3],res_16[[1]][3],res_18[[1]][3],
                               res_20[[1]][3],res_40[[1]][3],res_60[[1]][3])
  
  
  dataTwoTreeOne[[i]] <- rbind(res_06[[2]][3],res_08[[2]][3],
                               res_10[[2]][3],res_12[[2]][3],
                               res_14[[2]][3],res_16[[2]][3],res_18[[2]][3],
                               res_20[[2]][3],res_40[[2]][3],res_60[[2]][3])
  
  
  dataOneTreeTwo[[i]] <- rbind(res_06[[3]][3],res_08[[3]][3],
                               res_10[[3]][3],res_12[[3]][3],
                               res_14[[3]][3],res_16[[3]][3],res_18[[3]][3],
                               res_20[[3]][3],res_40[[3]][3],res_60[[3]][3])
  
  dataTwoTreeTwo[[i]] <- rbind(res_06[[4]][3],res_08[[4]][3],
                               res_10[[4]][3],res_12[[4]][3],
                               res_14[[4]][3],res_16[[4]][3],res_18[[4]][3],
                               res_20[[4]][3],res_40[[4]][3],res_60[[4]][3])
}


#### Lowess plot with 95% CI, I used "loess" instead of "lowess"
require(reshape)
require(ggplot2)

shapeClean<-function(TreeOne,TreeTwo){
  treeOnePoint <- sapply(TreeOne, t)
  treeTwoPoint <- sapply(TreeTwo, t)
  rownames(treeOnePoint) <- c(0.6,0.8,1,1.2,1.4,1.6,1.8,2,4,6)
  rownames(treeTwoPoint) <- c(0.6,0.8,1,1.2,1.4,1.6,1.8,2,4,6)
  
  meltTreeOne <- melt(treeOnePoint)[,-2]
  meltTreeTwo <- melt(treeTwoPoint)[,-2]
  
  meltTreeOnePlot <- meltTreeOne[order(meltTreeOne$X1),]
  meltTreeTwoPlot <- meltTreeTwo[order(meltTreeTwo$X1),]
  
  return(list(meltTreeOnePlot,meltTreeTwoPlot))
}

#meltSetOri<- rbind(cbind(meltTreeOnePlot,Group = "GeneTwoTreeOne"),cbind(meltTreeTwoPlot,Group = "GeneOneTreeTwo"))

##########
meltTreeOnePlot <- shapeClean(dataTwoTreeOne,dataOneTreeTwo)[[1]]
meltTreeTwoPlot <- shapeClean(dataTwoTreeOne,dataOneTreeTwo)[[2]]

plx1<-predict(loess(meltTreeOnePlot$value~meltTreeOnePlot$X1,span=2/3,degree=1,family="symmetric",iterations=4,surface="direct"), se=T)
plx2<-predict(loess(meltTreeTwoPlot$value~meltTreeTwoPlot$X1,span=2/3,degree=1,family="symmetric",iterations=4,surface="direct"), se=T)
pDiff <- ggplot() + geom_point(data = meltTreeOnePlot, aes(x = X1, y = value,color="GeneTwoTreeOne")) +
  geom_line(aes(x = meltTreeOnePlot$X1, y = plx1$fit, color = "GeneTwoTreeOne")) + 
  geom_line(aes(x = meltTreeOnePlot$X1, y = plx1$fit - qnorm(0.975)*plx1$se, color = "GeneTwoTreeOne"),linetype = 2) + 
  geom_line(aes(x = meltTreeOnePlot$X1, y = plx1$fit + qnorm(0.975)*plx1$se,  color = "GeneTwoTreeOne"), linetype = 2) +
  geom_point(data = meltTreeTwoPlot, aes(x = X1, y = value,color="GeneOneTreeTwo")) +
  geom_line(aes(x = meltTreeTwoPlot$X1, y = plx2$fit,color="GeneOneTreeTwo")) + 
  geom_line(aes(x = meltTreeTwoPlot$X1, y = plx2$fit - qnorm(0.975)*plx2$se,color="GeneOneTreeTwo"), linetype = 2) + 
  geom_line(aes(x = meltTreeTwoPlot$X1, y = plx2$fit + qnorm(0.975)*plx2$se,color="GeneOneTreeTwo"), linetype = 2) +
  # rd14, rd12 y(1, 1.7), rd10 y(1, 1.6)
  ylim(1, 1.6) + 
  theme(plot.margin = unit(c(0,0,0,0), "lines"),plot.title = element_text(hjust = 0.5),
        plot.background = element_blank()) +
  ggtitle("Under Different Tree") +
  labs(x = "C Value", y = "Ratio", color = "") + theme(legend.position = "bottom",legend.text = element_text( size=12))
# Renew data
meltTreeOnePlotS <- shapeClean(dataOneTreeOne,dataTwoTreeTwo)[[1]]
meltTreeTwoPlotS <- shapeClean(dataOneTreeOne,dataTwoTreeTwo)[[2]]

plx1S<-predict(loess(meltTreeOnePlotS$value~meltTreeOnePlotS$X1,span=2/3,degree=1,family="symmetric",iterations=4,surface="direct"), se=T)
plx2S<-predict(loess(meltTreeTwoPlotS$value~meltTreeTwoPlotS$X1,span=2/3,degree=1,family="symmetric",iterations=4,surface="direct"), se=T)
pSame <- ggplot() + geom_point(data = meltTreeOnePlotS, aes(x = X1, y = value,color="GeneOneTreeOne")) +
  geom_line(aes(x = meltTreeOnePlotS$X1, y = plx1S$fit, color = "GeneOneTreeOne")) + 
  geom_line(aes(x = meltTreeOnePlotS$X1, y = plx1S$fit - qnorm(0.975)*plx1S$se, color = "GeneOneTreeOne"),linetype = 2) + 
  geom_line(aes(x = meltTreeOnePlotS$X1, y = plx1S$fit + qnorm(0.975)*plx1S$se,  color = "GeneOneTreeOne"), linetype = 2) +
  geom_point(data = meltTreeTwoPlotS, aes(x = X1, y = value,color="GeneTwoTreeTwo")) +
  geom_line(aes(x = meltTreeTwoPlotS$X1, y = plx2S$fit,color="GeneTwoTreeTwo")) + 
  geom_line(aes(x = meltTreeTwoPlotS$X1, y = plx2S$fit - qnorm(0.975)*plx2S$se,color="GeneTwoTreeTwo"), linetype = 2) + 
  geom_line(aes(x = meltTreeTwoPlotS$X1, y = plx2S$fit + qnorm(0.975)*plx2S$se,color="GeneTwoTreeTwo"), linetype = 2) +
  # rd14, rd12 y(1, 1.7), rd10 y(1, 1.6)
  ylim(1, 1.6) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "lines"),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank()) +
  ggtitle("Under Same Tree")+
  labs(x = "C Value", y = "Ratio", color = "") + theme(legend.position = "bottom",legend.text = element_text( size=12))

library(gridExtra)
grid.arrange(pDiff,pSame,ncol=2)
