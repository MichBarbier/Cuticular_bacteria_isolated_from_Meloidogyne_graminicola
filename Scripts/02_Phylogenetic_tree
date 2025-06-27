################################################
# This script allows to reproduce the Figure 2 #

setwd(C:/[[YOUR PATH]]/Cuticular_bacteria_isolated_from_Meloidogyne_graminicola)

######## Libraries
library(ape)
library(Biostrings)
library(ggtree)
library(msa)
library(tidyverse)

######## Alignment 
## Load sequence corrected with BioEdit software as ".fasta" and alignin it
Sequences <- readDNAStringSet(file = "Data/Sequences_for_tree.fas")
Sequences[1:17] <- reverseComplement(Sequences[1:17])

## Alignment using the Muscle algorithm
Alignment <- msaMuscle(inputSeqs = Sequences, cluster = "neighborjoining", type = "dna", maxiters = 100, verbose = TRUE)
Alignment_matrix <- as.matrix(Alignment)
Alignment_DNAbin <- as.DNAbin(Alignment_matrix)

######## Phylogenetic tree
## Create a tree using the Neighbor joining algorithm and the Kimura 80 model
Distances <- dist.dna(Alignment_DNAbin, model = "K80")
Tree <- nj(Distances)

## Use the sequence "NR_114042.1 Escherichia coli" as root
Tree <- root(Tree, outgroup = "NR_114042.1 Escherichia coli", r = TRUE)

## 10000 Bootstraps 
Boots <- boot.phylo(Tree, Alignment_DNAbin, FUN = function(Tree) nj(dist.dna(Tree)), B = 10000)

## Allow to know the number of tip and node to plot the bootstrap scores 
Ntip(Tree)+1
Ntip(Tree)+Nnode(Tree)

## Create a data frame with the boostrap scores, the row numbers represent the nodes where each scores will be ploted
## "61:119" are respectively the results of "Ntip(Tree)+1" and "Ntip(Tree)+Nnode(Tree)"
Boots_scores <- data.frame()
Boots_scores[61:119,"Boots"] <- as.data.frame(round(Boots/10000*100)) # Bootstrap scores are transformed as percentage

## Plot the phylogenetic tree
p1 <- ggtree(Tree, layout = "rectangular", options(ignore.negative.edge = TRUE))+
  geom_tiplab(size = 2)+
  geom_nodepoint(color = "black", size = 0.5)+
  geom_nodelab(aes(label = Boots_scores$Boots), nudge_x = -0.002, size = 1.75, color = "red", fontface = "bold")+
  geom_treescale(y = -3, fontsize = 3)
p1

## The Figure 2 was modified using Microsoft PowerPoint
## Colors were added depending of the isolation soil and overlaps of the bootstrap scores and the tree were removed 
ggsave(plot = p1, dpi = 1000, device = "svg", width = 12, height = 15, filename = "Results/Figure_2.svg")

