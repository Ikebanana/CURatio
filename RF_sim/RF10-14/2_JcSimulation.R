#####
#####
#####

library(ape)
library(phangorn)

#setwd('//as-phoenix1.ad.uky.edu/Statistics/qka222/Desktop/r/08052016RFDistanceAndSim/rd10_dna/')
setwd("E:/work/r/08052016RFDistanceAndSim/rd10_dna/")
cValue<-c("06","07","08","09","10","12","14","16","18","20","40","60")
sNum<-c(0:9)
for(i in cValue){
  for(j in sNum){
    ## Read the data set and set work station
    fileDir<-paste("depth",i,"/set",j,sep="")
    fileList<-list.files(path=fileDir,pattern = ".phylip")
    treeDir<-paste("data/depth",i,"/set",j,"/Sp",i,".nex",sep = "")
    outDir<-paste("data/depth",i,"/set",j,"/",sep = "")
    # Compared to JcSimulation3.0.R, the only different part is here
    # I change the code from "read.nexus" to "read.tree" 
    Tree<-read.tree(treeDir)
    spTreeOne<-compute.brlen(Tree[[1]])
    spTreeTwo<-compute.brlen(Tree[[2]])
    
    ## Begin to calculate the ratio 
    geneOneSpOne<-simJC(spTreeOne,fileList[1:1000],i,j)
    geneOneSpTwo<-simJC(spTreeTwo,fileList[1:1000],i,j)
    geneTwoSpOne<-simJC(spTreeOne,fileList[1001:2000],i,j)
    geneTwoSpTwo<-simJC(spTreeTwo,fileList[1001:2000],i,j)
    
#     RF.dist(spTreeOne,spTreeTwo)
#     plot(spTreeOne)
#     plot(spTreeTwo)
#     ks.test(geneOneSpOne,geneTwoSpOne)
#     ks.test(geneOneSpTwo,geneTwoSpTwo)
    ## Write the table we want
    # Gene One Tree One
    result.mix11<-data.frame(fileList[1:1000],geneOneSpOne,stringsAsFactors=FALSE)
    colnames(result.mix11)<-c('Names','Ratio')  
    write.table(result.mix11,paste(outDir,'All_1.0_tree1.txt',sep=""),sep='\t')
    # Gene One Tree Two
    result.mix12<-data.frame(fileList[1:1000],geneOneSpTwo,stringsAsFactors=FALSE)
    colnames(result.mix12)<-c('Names','Ratio')  
    write.table(result.mix12,paste(outDir,'All_1.0_tree2.txt',sep=""),sep='\t')
    # Gene Two Tree One
    result.mix21<-data.frame(fileList[1001:2000],geneTwoSpOne,stringsAsFactors=FALSE)
    colnames(result.mix21)<-c('Names','Ratio')  
    write.table(result.mix21,paste(outDir,'All_2.0_tree1.txt',sep=""),sep='\t')
    # Gene Two Tree Two
    result.mix22<-data.frame(fileList[1001:2000],geneTwoSpTwo,stringsAsFactors=FALSE)
    colnames(result.mix22)<-c('Names','Ratio')  
    write.table(result.mix22,paste(outDir,'All_2.0_tree2.txt',sep=""),sep='\t')
    
    ## Histgrom
    histGG(i,j)
  }
}

simJC<-function(spTree,fileList,cValue,sNum){
  #### The sum of branch length of each tree 
  b<-c()
  B<-c()
  output<-c()
  stan.tree<-spTree

  #### Calculating the ratio
  for(i in 1:1000){
    workDir<-paste("depth",cValue,"/set",sNum,"/",fileList[[i]],sep="")
    # STEP 1 
    # Reading the original tree, we also need to clean the name value, it 
    # is kind of messy
    data<-read.dna(workDir)
    data.list<-phyDat(data,type='DNA')
    
    # STEP 2 
    # Reading the MLE tree
    dm<-dist.hamming(data.list)
    treeNJ<-NJ(dm)
    fit<- pml(treeNJ,data.list)
    treeJC <- optim.pml(fit, optNni=TRUE, model="JC")
    b[i]<-sum(treeJC$tree$edge.length)
    
    # STEP 3
    
    fit2<- pml(stan.tree,data.list)
    fit.opt<- optim.pml(fit2, optNni=FALSE, model="JC") ### MLE tree under the JC model with constraint
    B[i]<-sum(fit.opt$tree$edge.length)
    #OUTPUT
    output[i]<-B[i]/b[i]
  }
  return(output)
}

histGG<-function(cValue,sNum){
  require(ggplot2)
  workDir<-paste("data/depth",cValue,"/set",sNum,"/",sep = "")
  
  fileList<-list.files(workDir,'.')#
  geneOneTreeOne<-as.data.frame(read.table(paste(workDir,fileList[1],sep=""))[,2])
  geneOneTreeTwo<-as.data.frame(read.table(paste(workDir,fileList[2],sep=""))[,2])
  geneTwoTreeOne<-as.data.frame(read.table(paste(workDir,fileList[3],sep=""))[,2])
  geneTwoTreeTwo<-as.data.frame(read.table(paste(workDir,fileList[4],sep=""))[,2])
  
  geneOneTreeOne$Type<-'GeneOne'
  geneOneTreeTwo$Type<-'GeneOne'
  geneTwoTreeOne$Type<-'GeneTwo'
  geneTwoTreeTwo$Type<-'GeneTwo'
  
  colnames(geneOneTreeOne)[1]<-'Ratio'
  colnames(geneOneTreeTwo)[1]<-'Ratio'
  colnames(geneTwoTreeOne)[1]<-'Ratio'
  colnames(geneTwoTreeTwo)[1]<-'Ratio'
  ####Under Tree One
  
  treeOne <- rbind(geneOneTreeOne,geneTwoTreeOne)
  treeTwo <- rbind(geneOneTreeTwo,geneTwoTreeTwo)
  
  
  jpeg(filename = paste(workDir,'TreeOne',cValue,'.png',sep=""), width = 400, height = 300)
  picOne<-ggplot(treeOne,  aes(Ratio, fill = Type))+ 
    geom_histogram(data = subset(treeOne, Type == 'GeneOne'),  alpha = 0.2,  binwidth = 0.05)+ 
    geom_histogram(data = subset(treeOne, Type == 'GeneTwo'),  alpha = 0.2,  binwidth = 0.05)+
    scale_fill_manual(values = c("red", "blue"))+
    geom_vline(xintercept = quantile(subset(treeOne, Type == 'GeneOne')[,1],.95), colour = 'red', alpha =0.4)+
    geom_vline(xintercept = quantile(subset(treeOne, Type == 'GeneTwo')[,1],.95), colour = 'blue', alpha =0.4)+
    xlim(c(1,3))+
    ggtitle('Histogram of Ratios\n Tree One')
  print(picOne)
  dev.off()
  
  jpeg(filename = paste(workDir,'TreeTwo',cValue,'.png',sep=""), width = 400, height = 300)
  picTwo<-ggplot(treeOne,  aes(Ratio, fill = Type))+ 
    geom_histogram(data = subset(treeTwo, Type == 'GeneOne'),  alpha = 0.2, binwidth = 0.05)+ 
    geom_histogram(data = subset(treeTwo, Type == 'GeneTwo'),  alpha = 0.2,  binwidth = 0.05)+
    scale_fill_manual(values = c("red", "blue"))+
    geom_vline(xintercept = quantile(subset(treeTwo, Type == 'GeneOne')[,1],.95), colour = 'red', alpha =0.4)+
    geom_vline(xintercept = quantile(subset(treeTwo, Type == 'GeneTwo')[,1],.95), colour = 'blue', alpha =0.4)+
    xlim(c(1,3))+
    ggtitle('Histogram of Ratios\n Tree Two')
  print(picTwo)
  dev.off()
}
