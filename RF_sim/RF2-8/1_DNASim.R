##### Small RF disntance
##### Qiwen Kang
##### 01/08/2018


# Description -------------------------------------------------------------
# This code is for creating the DNA alignments based on the gene
# trees we already have. Since we already include the DNA alignment
# in our data sets, you don't need this program any more.

# Setting up --------------------------------------------------------------
library(ape)
library(phangorn)
setwd("/home/qiwen/r/180108SmallRFDistance/")

path <- "./data/rd8/depth06/set0/"

# Function ----------------------------------------------------------------
printf <- function(...) invisible(print(sprintf(...)))

mdk <- function(path){
  fileNames <- list.files(path)
  gene1 <- read.nexus(paste(path,fileNames[1],sep=""))
  gene2 <- read.nexus(paste(path,fileNames[2],sep=""))
  for(i in 1:1000){
    gene1[[i]]$edge.length <- gene1[[i]]$edge.length/sum(gene1[[i]]$edge.length)
    gene2[[i]]$edge.length <- gene2[[i]]$edge.length/sum(gene2[[i]]$edge.length)
    
    wt1 <- printf(paste(path, "GeneOne%d.tre", sep = ""), i)
    write.tree(gene1[[i]], wt1, append = FALSE, digits = 4)
    comd <- printf("./run_paml_JC %s > jnk2", wt1)
    system(comd)
    
    wt2 <- printf(paste(path, "GeneTwo%d.tre", sep = ""), i)
    write.tree(gene2[[i]], wt2, append = FALSE, digits = 4)
    comd <- printf("./run_paml_JC %s > jnk2", wt2)
    system(comd)
  }
}

# Main --------------------------------------------------------------------
disName <- c("rd2","rd4","rd6","rd8")
depth <- paste("depth",c("06","08","10","12","14",
                         "16","18","20","40","60"),sep="")
setName <-paste("set",0:9,sep="") 

dirFrame <- expand.grid(disName, depth, setName)

workDir <- paste("./data",dirFrame$Var1,dirFrame$Var2,dirFrame$Var3,"",sep="/")

lapply(workDir, mdk)

# for(i in 1:400) mdk(workDir[i])
