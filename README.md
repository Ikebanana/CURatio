# CURatio

## Description
All the programs in this repository is for CURatio, which uses ratios of total branch lengths in gene trees to help identify phylogenetic outliers in a given set of ortholog groups from multiple genomes. An advantage of CURatio over other methods is that genes absent from and/or duplicated in some genomes can be included in the analysis. The constraint tree is the 65% consensus of presumed housekeeping genes (those with orthologues in all 12 genomes) on the JC model. For more details, please check *CURatio: Genome-wide phylogenomic analysis method using ratios of total branch lengths*.

## Algorithm
1. Calculating the sum of branch lengths of the tree without consensus tree, **b**.
  - Performing the neighbor-joining tree estimation of Saitou and Nei (1987).
  - Computing the likelihood of a phylogenetic tree given a sequence alignment and JC model.
  - Optimizing the different model parameters. We could change the shape of the phylogenentic tree.

2. Calculating the sum of branch lengths of the tree with consensus tree, **B**.
  - Inputting the consensus tree as the phylogenentic tree.
  - Computing the likelihood of the phylogenetic tree given a sequence alignment and JC model.
  - Optimizing the different model parameters. Do **NOT** change the shape of the phylogenentic tree.

3. Calculating the ratio 

      **Ratio = B/b**
      
## R package required for CURatio
* ape
* phangorn

## Example
```{r}
library(ape)
library(phangorn)

phyData <- read.phyDat(file="./data/eas-ccl-wg_10006_2_2_1_2_2-dna-trimmed.fasta",format='fasta',type='DNA')
consTree <- read.tree("./consensusTree_0.65.txt")

ratio <- CURatio(phyData, consTree)
```

## Relevant Citations
Saitou, Naruya, and Masatoshi Nei. "The neighbor-joining method: a new method for reconstructing phylogenetic trees." Molecular biology and evolution 4.4 (1987): 406-425.
