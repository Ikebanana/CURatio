##### Outlier group DNA generating
##### Qiwen Kang
##### 01/16/2018

# Functions ----------------------------------------------------------------
# The first two functions are for simulating the DNA aligments
printf <- function(...) invisible(print(sprintf(...)))
mdk <- function(path){
  fileNames <- list.files(path)
  gene1 <- read.nexus(paste(path,fileNames[1],sep=""))
  
  for(i in 1:6000){
    gene1[[i]]$edge.length <- gene1[[i]]$edge.length/sum(gene1[[i]]$edge.length)
    
    wt1 <- printf(paste(path, "outlier%d.tre", sep = ""), i)
    write.tree(gene1[[i]], wt1, append = FALSE, digits = 4)
    comd <- printf("./run_paml_JC %s > jnk2", wt1)
    system(comd)
  }
}

# These two functions are for calculating CURatio
sevenJC<-function(spTree,fileDir){
  #### The sum of branch length of each tree 
  b <- rep(NA,6000)
  B <- rep(NA,6000)
  output <- rep(NA,6000)
  stan.tree <- compute.brlen(spTree)
  fileList <- list.files(path=fileDir,pattern = ".phylip")
  
  #### Calculating the ratio
  for(i in 1:6000){
    workDir<-paste(fileDir,fileList[[i]],sep="")
    # STEP 1 
    # Reading the original tree, we also need to clean the name value, it 
    # is kind of messy
    dataOri<-read.dna(workDir)
    data <- as.list(dataOri)
    data$aa <- as.list(dataOri)$a
    data$ee <- as.list(dataOri)$e
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

sevenJC_out<-function(spTree,fileDir){
        #### The sum of branch length of each tree 
        b <- rep(NA,6000)
        B <- rep(NA,6000)
        output <- rep(NA,6000)
        stan.tree <- compute.brlen(spTree)
        fileList <- list.files(path=fileDir,pattern = ".phylip")
        
        #### Calculating the ratio
        for(i in 1:6000){
                workDir<-paste(fileDir,fileList[[i]],sep="")
                # STEP 1 
                # Reading the original tree, we also need to clean the name value, it 
                # is kind of messy
                dataOri<-read.dna(workDir)
                data <- as.list(dataOri)
                data$aa <- as.list(dataOri)$e
                data$ee <- as.list(dataOri)$e
                data$e <- as.list(dataOri)$a
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


# Setting up --------------------------------------------------------------
library(ape)
library(phangorn)
setwd("~/r/180116sevenLeaves/")

# Creating DNA alignment (Just running 1 time) --------------------------------------------------------------------
path_non <- "./data2/non_outlier/"
path_out <- "./data2/out/"
# mdk(path_non)
# mdk(path_out)

# Calculating CURatio -----------------------------------------------------
species <- unroot(read.nexus("./data2/sp30.nex")[[30]])

# Non_outlier with in-paralogs
output_non <- sevenJC(species,path_non)

fileList <- list.files(path=path_non,pattern = ".phylip")
result.mix1<-data.frame(fileList,output_non,stringsAsFactors=FALSE)
colnames(result.mix1)<-c('Names','Ratio')  
write.table(result.mix1,paste(path_non,'Ratio_non.txt',sep=""),sep='\t')

# out_paralogs
output_out <- sevenJC_out(species,path_non)

fileList2 <- list.files(path=path_non,pattern = ".phylip")
result.mix2<-data.frame(fileList2,output_out,stringsAsFactors=FALSE)
colnames(result.mix2)<-c('Names','Ratio')  
write.table(result.mix2,paste(path_non,'Ratio_out.txt',sep=""),sep='\t')

non <- read.table("./data2/non_outlier/Ratio_non.txt")
out <- read.table("./data2/non_outlier/Ratio_out.txt")
non[,"Group"] <- "Non_outlier"
out[,"Group"] <- "Out_paralog"
treeOne <- rbind(non,out)
require(ggplot2)    
jpeg(filename = paste(workDir,'TreeOne',cValue,'.png',sep=""), width = 400, height = 300)
picOne <- ggplot(treeOne,  aes(Ratio, fill = Group))+ 
          geom_histogram(data = subset(treeOne, Group == 'Non_outlier'),  alpha = 0.2,  binwidth = 0.05)+ 
          geom_histogram(data = subset(treeOne, Group == 'Out_paralog'),  alpha = 0.2,  binwidth = 0.05)+
          theme(plot.title = element_text(hjust = 0.5))+
          scale_fill_manual(values = c("red", "blue"))+
          xlim(c(0.9,2.5))+
          ylab("Count")+
          ggtitle('Histogram of Ratios')
print(picOne)
dev.off()

