# Description -------------------------------------------------------------
# In this script, we just include our CURatio functions
# The first one "CURatio" is for calculating CURatios.
# The second one is for the case when we have multi-representatives
# of one gene. The "dupl" function could duplicate the a spacific gene
# in a given species tree.

CURatio<-function(stdTree){
  library(ape)
  library(phangorn)
  # In this part, we define some variables to save the value
  b<-c() # It is the sum of branch length of the tree without a consensus tree.
  # You can change the topology of the tree.
  B<-c() # It is the sum of branch length of the tree with a consensus tree.
  # You cannot change the topology of the tree.
  dist<-c() # It is the Robinson-Foulds distance
  output<-c() # This is the ratios
  fileListNew<-c() # This is the new name list to save the name we need.
  # Since we want to remove the 
  # Given a consensus tree, this stdTree is getted from the user. 
  # Since the tree have no branch lengths, so we calculate the branch lengths
  # using Grafen's method (Grafen, 1989).
  stan.tree<-compute.brlen(stdTree) 
  # The tip labels of the consensus tree
  stan.name<-Tree$tip.label
  # Reading data set from the current work station
  fileList<-list.files(path='.',pattern='.fasta')
  
  # From here, we begin our for loop, reading each DNA sequence into RAM 
  # and calculating the ratios. 
  for(i in 1:length(fileList)){
    # Reading the DNA sequence into RAM
    data.list<-read.phyDat(file=fileList[[i]],format='fasta',type='DNA')
    # Since we want to change the name of the tip labels
    splitValue<-sapply(names(data.list)[1:length(data.list)],function(x) strsplit(x,"|",fixed=T))
    nameValue<-lapply(splitValue,function(x) x[1])
    names(data.list)<-unlist(nameValue)
    
    # We separate the DNA sequence into 3 different cases: <5, 5~11, 12.
    if (5<=length(data.list)&& length(data.list)<12){
      
      # STEP 1: Calculating the sum of branch lengths of the tree without consensus tree
      # Computing pairwise distances for an object of class phyDat.
      dm<-dist.hamming(data.list)
      # Performing the neighbor-joining tree estimation of Saitou and Nei (1987).
      treeNJ<-NJ(dm)
      # Computing the likelihood of a phylogenetic tree given a sequence alignment and a model.
      fit<- pml(treeNJ,data.list)
      # Optimizing the different model parameters.
      treeFit <- optim.pml(fit, optNni=TRUE, model="JC")
      # 
      b[i]<-sum(treeFit$tree$edge.length)
      
      # STEP 2: Calculating the sum of branch lengths of the tree with consensus tree
      # If the number of tips is less than 12, we need to remove some of the missing
      # tips from the consensus tree.
      data.name<-attr(data.list,'names')
      index<- which(stan.name%in%data.name)
      remove.name<-stan.name[-index]
      tree.new<-drop.tip(stan.tree,remove.name)
      
      # Computing the likelihood of a phylogenetic tree given a sequence alignment and a model.
      fit2<- pml(tree.new,data.list)
      # Optimizing the different model parameters.
      treeFit2<- optim.pml(fit2, optNni=FALSE, model="JC") ### MLE tree under the JC model with constraint
      # 
      B[i]<-sum(treeFit2$tree$edge.length)
      
      # We calculate the Robinson-Foulds distance to compare the two trees
      dist[i]<-RF.dist(treeFit$tree,treeFit2$tree)
      # We keep all the new names of the DNA sequences
      fileListNew[i]<-paste(unlist(strsplit(fileList[i],"-wg"))[1],unlist(strsplit(fileList[i],"-wg"))[2],sep='')
      # Calculating the ratios
      output[i]<-B[i]/b[i]
      
    } else if(length(data.list)==12){
      # Here, we do the same work. But the case is when the number of
      # tree's tips is equal to 12.
      # STEP 1 
      dm<-dist.hamming(data.list)
      treeNJ<-NJ(dm)
      fit<- pml(treeNJ,data.list)
      treeFit <- optim.pml(fit, optNni=TRUE, model="JC")
      b[i]<-sum(treeFit$tree$edge.length)
      
      # STEP 2
      
      fit2<- pml(stan.tree,data.list)
      treeFit2<- optim.pml(fit2, optNni=FALSE, model="JC") ### MLE tree under the JC model with constraint
      B[i]<-sum(treeFit2$tree$edge.length)
      
      dist[i]<-RF.dist(treeFit$tree,treeFit2$tree)
      fileListNew[i]<-paste(unlist(strsplit(fileList[i],"-wg"))[1],unlist(strsplit(fileList[i],"-wg"))[2],sep='')
      #OUTPUT
      output[i]<-B[i]/b[i]
    } else next
    
  }
  # Since the numbers of some of the trees's tips are less than 5, so we get some
  # NA value in our data set. We use "position" to mark the position.
  position<-is.na(output)
  # We remove the NA value from the data set.
  result.mix<-data.frame(fileListNew[!position],output[!position],dist[!position],stringsAsFactors=FALSE)
  # We return the final dataframe data set at last.
  return(result.mix)
}

dupl <- function(tip_label,consen_path){
  # tip_label: The tip lable you want to duplicate
  # consen_path: The directory of your consensus tree.
  if(require(ape)){
    dupl_text <- readLines(consen_path)
    if(grepl(tip_label, dupl_text)){
      new_pattern <- paste("(",tip_label,",",tip_label,")",sep = "")
      new_lines <- gsub(tip_label,new_pattern,dupl_text)
      conTree <- read.tree(text = new_lines)
      return(conTree)
    }
    else{
      warning("No tip label matched in the tree.")
    }
  }
  else{
    warning("This function requires 'ape' package.")
  }
}
