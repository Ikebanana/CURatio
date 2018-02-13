##### Robinson-Foulds distance and Simulation
##### Qiwen Kang
##### 08/05/2016

#### Details
### This code is for calculating the distance between two trees.
### I got 2000 trees in Mesquite and tried to class these trees into 3 different groups

library(ape)
library(phangorn)
mainDir<-"//as-phoenix1.ad.uky.edu/Statistics/qka222/Desktop/r/08052016RFDistanceAndSim"
setwd(mainDir)
#### Creating the Dirctories
subDir<-c("set0","set1","set2","set3","set4","set5","set6","set7","set8","set9")
cValue<-c("06","07","08","09","10","12","14","16","18","20","40","60")
newDir<-paste("./data/rd10/depth",cValue,sep="")
for(i in 1:12){
  dir.create(newDir[i])
}

for(i in cValue){
  newDir<-paste("./data/rd10/depth",i,"/",subDir,sep="")
  for(j in 1:10){
    dir.create(newDir[j])
  }
}



#### Read data set and find the specific distance
spTree<-function(cValue){
  
  tree<-read.nexus(paste("spe",cValue,".nex",sep=""))
  rd<-matrix(,2000,2000)
  
  for(i in 1:1999){
    for(j in (i+1):2000){
      rd[i,j]<-RF.dist(tree[[i]],tree[[j]])
    }
  }
  
  ### pick up the distance you want.
  
  rd10<-which(rd==10)
  rd12<-which(rd==12)
  rd14<-which(rd==14)
  
  #### Output the distance
  ### Distance 10
  row10<-ceiling(rd10%%2000)
  coloum10<-ceiling(rd10/2000)
  ranNum10<-sample(1:length(coloum10),10)
  for(i in 1:10){
    write.tree(tree[c(row10[ranNum10[i]],coloum10[ranNum10[i]])],file = paste("./data/rd10/depth",cValue,"/set",i-1,"/Sp",cValue,".nex",sep=""))
    
  }

  
  ### Distance 12
  row12<-ceiling(rd12%%2000)
  coloum12<-ceiling(rd12/2000)
  ranNum12<-sample(1:length(coloum12),10)
  for(i in 1:10){
    write.tree(tree[c(row12[ranNum12[i]],coloum12[ranNum12[i]])],file = paste("./data/rd12/depth",cValue,"/set",i-1,"/Sp",cValue,".nex",sep=""))
  }

  
  ### Distance 14
  row14<-ceiling(rd14%%2000)
  coloum14<-ceiling(rd14/2000)
  ranNum14<-sample(1:length(coloum14),10)
  for(i in 1:10){
    write.tree(tree[c(row14[ranNum14[i]],coloum14[ranNum14[i]])],file = paste("./data/rd14/depth",cValue,"/set",i-1,"/Sp",cValue,".nex",sep=""))
  }
  
  colTable<-rbind(row10[ranNum10],coloum10[ranNum10],row12[ranNum12],
        coloum12[ranNum12],row14[ranNum14],coloum14[ranNum14])
  row.names(colTable)<-c("r10","c10","r12","c12","r14","c14")
  write.table(colTable,paste("ranNum",cValue,".txt",sep=""),sep="\t")
  return(colTable)
}



