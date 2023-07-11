install.packages('ape')
install.packages('phangorn')
install.packages('seqinr')
install.packages('phylotools')
install.packages("TreeTools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("seqLogo")

library(ape)
library(phangorn)
library(seqinr)
library(phylotools)
library(TreeTools)

#clear existing environment
rm(list = ls())

#set the directory - alter according to your file structure
setwd("")

# read the newick file and the FASTA sequences
tree <-  ape::read.tree('final-project/data_files/tree_file.nwk')
phySeq <- read.phyDat('final-project/data_files/HA_cds.fa', format = 'fasta', type = "DNA")

# creates a phylo object
anc.acctran <- acctran(tree, phySeq)
anc.acctran.othermethod.new  <- di2multi(anc.acctran)

# creates phydat object --> replace tree with the phylo object to be used
anc.acctran.othermethod.2 <- ancestral.pars(anc.acctran.othermethod.new, phySeq, "ACCTRAN")

# save phylo object
output_fitch_tree_othermethod <- as.Newick(anc.acctran.othermethod.new)
ape::write.tree(anc.acctran.othermethod.new, file='output_tree_reconstructed.txt')
ape::write.nexus(anc.acctran.othermethod.new, file='output_tree_reconstructed.nex')

# save info from phylo object 
edge_connections <- data.frame(anc.acctran.othermethod.new$edge)
write.csv(edge_connections, "output_tree_reconstructed_phylo_edge.csv", row.names=FALSE, quote=FALSE) 
edge_length <- data.frame(anc.acctran.othermethod.new$edge.length)
write.csv(edge_length, "output_tree_reconstructed_phylo_edge-length.csv", row.names=FALSE, quote=FALSE) 
node_label <- data.frame(anc.acctran.othermethod.new$node.label)
write.csv(node_label, "output_tree_reconstructed_phylo_node-label.csv", row.names=FALSE, quote=FALSE) 
tip_label <- data.frame(anc.acctran.othermethod.new$tip.label)
write.csv(tip_label, "output_tree_reconstructed_phylo_tip-label.csv", row.names=FALSE, quote=FALSE) 
nnode <- data.frame(anc.acctran.othermethod.new$Nnode)
write.csv(nnode, "output_tree_reconstructed_phylo_nnode.csv", row.names=FALSE, quote=FALSE) 

# plot tree with reconstructed sequences
plotAnc(phyTree, anc.acctran, 100)
title("ACCTRAN")

# parse the phyDat object to obtain sequences and corresponding node name/number
# set empty data structures required
sequences <- rep(c(''), times = length(anc.acctran.othermethod.2[,2]))
anc.acctran.othermethod.df <- as.data.frame(anc.acctran.othermethod.2)
seq_names <- colnames(anc.acctran.othermethod.df)

phyDat_df <- data.frame(seq_names, sequences)

# loop over the first 1625 sequences (these are the leaves of tree)
for (seqnum in 1:1625){
  # obtain sequence matrix
  seq <- anc.acctran.othermethod.2[[seq_names[seqnum]]]
  seq <- as.data.frame((seq))
  # replace numbers with relevant letters
  seq$a[seq$a > 0] <- 'a'
  seq$c[seq$c > 0] <- 'c'
  seq$g[seq$g > 0] <- 'g'
  seq$t[seq$t > 0] <- 't'

  sequence <- '' 
  
  # use this updated matrix to reconstruct sequence at that node 
  for (i in 1:(length(seq$a))){
    row <- as.list((seq[i,]))
    boolrow <- which(row != '0', arr.ind = T)
    letter <- sample(row[boolrow], 1)
    sequence <- paste(sequence, letter)
  }
  
  # remove spaces from sequence
  sequence <- gsub(" ","",sequence) 
  
  # append sequence to data frame
  phyDat_df$sequences[seqnum] <- sequence
}

# repeat for rest of sequences (corresponding to intermediate nodes)
for (seqnum in 1626:length(seq_names)){
  seq <- anc.acctran.othermethod.2[[seq_names[seqnum]]]
  seq <- as.data.frame((seq))
  seq$V1[seq$V1 > 0] <- 'a'
  seq$V2[seq$V2 > 0] <- 'c'
  seq$V3[seq$V3 > 0] <- 'g'
  seq$V4[seq$V4 > 0] <- 't'
  
  sequence <- '' 
  
  for (i in 1:(length(seq$V1))){
    row <- as.list((seq[i,]))
    boolrow <- which((row != '0'), arr.ind = T)
    letter <- sample(row[boolrow], 1)
    sequence <- paste(sequence, letter)
  }
  
  sequence <- gsub(" ","",sequence) 
  
  phyDat_df$sequences[seqnum] <- sequence
}

# save the data frame as a CSV for future use
write.csv(phyDat_df, "reconstructed_sequences_phyDat_ramdomly-sampled.csv", row.names=FALSE, quote=FALSE) 

  
