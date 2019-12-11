  
CURatio <- function(phyData, consTree){
  stan.tree <- compute.brlen(consTree) 
  # The tip labels of the consensus tree
  stan.name <- consTree$tip.label
  
  if(length(phyData) > length(stdTree$tip.label)){
    ratio <- NA
  }
  else if(length(phyData) == length(stdTree$tip.label)){
    dm <- dist.hamming(phyData)
    treeNJ <- NJ(dm)
    fit <- pml(treeNJ,phyData)
    treeFit <- optim.pml(fit, optNni=TRUE, model="JC")
    b <- sum(treeFit$tree$edge.length)
    
    # STEP 2
    
    fit2 <- pml(stan.tree,phyData)
    treeFit2 <- optim.pml(fit2, optNni=FALSE, model="JC") ### MLE tree under the JC model with constraint
    B <- sum(treeFit2$tree$edge.length)
    
    ratio <- B/b
  }
  else if(2 < length(phyData) && length(phyData) < length(stdTree$tip.label)){
    dm <- dist.hamming(phyData) 
    treeNJ <- NJ(dm)
    # Computing the likelihood of a phylogenetic tree given a sequence alignment and a model.
    fit <- pml(treeNJ,phyData)
    # Optimizing the different model parameters.
    treeFit <- optim.pml(fit, optNni=TRUE, model="JC")
    # 
    b <- sum(treeFit$tree$edge.length)
    
    data.name <- attr(phyData,'names')
    index <- which(stan.name %in% data.name)
    remove.name <- stan.name[-index]
    tree.new <- drop.tip(stan.tree,remove.name)
    
    # Computing the likelihood of a phylogenetic tree given a sequence alignment and a model.
    fit2 <- pml(tree.new,phyData)
    # Optimizing the different model parameters.
    treeFit2 <- optim.pml(fit2, optNni=FALSE, model="JC") ### MLE tree under the JC model with constraint
    # 
    B<-sum(treeFit2$tree$edge.length)
    
    ratio = B/b
  }
  else {ratio <- NA}
  
  return(ratio)
}