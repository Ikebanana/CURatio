# CURatio

## Description
These are the coded we used to calculate the ratios. The constraint tree we used is the 65% consensus of presumed housekeeping genes (those with orthologues in all 12 genomes) on the JC model.

## Algorithm
1. Calculating the sum of branch lengths of the tree without consensus tree, **b**.
  - Performing the neighbor-joining tree estimation of Saitou and Nei (1987).
  - Computing the likelihood of a phylogenetic tree given a sequence alignment and JC model.
  - Optimizing the different model parameters. 

2. Calculating the sum of branch lengths of the tree with consensus tree, **B**.
  - Inputting the consensus tree as the phylogenentic tree.
  - Computing the likelihood of the phylogenetic tree given a sequence alignment and JC model.
  - Optimizing the different model parameters. 

3. Calculating the ratio
  $$Ratio = \frac{B}{b}$$
