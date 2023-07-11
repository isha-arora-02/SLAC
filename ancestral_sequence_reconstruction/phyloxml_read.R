install.packages('ape')
install.packages('phylotools')
install.packages("TreeTools")
install.packages("devtools")
devtools::install_github("USCBiostats/rphyloxml")

library(rphyloxml)
library(ape)

#clear existing environment
rm(list = ls())

#set the directory
setwd("{add rest of path here and remove brackets}/final-project-main/data_files")

# read xml file
xmltree <- read_phyloxml('phy_tree.xml')

phylo_obj <- phyloxml2phylo(xmltree)

# save info from phylo object 
edge_connections <- data.frame(phylo_obj[[1]]$edge)
write.csv(edge_connections, "output_tree_reconstructed_phylo_edge.csv", row.names=FALSE, quote=FALSE) 
edge_length <- data.frame(phylo_obj[[1]]$edge.length)
write.csv(edge_length, "output_tree_reconstructed_phylo_edge-length.csv", row.names=FALSE, quote=FALSE) 
tip_label <- data.frame(phylo_obj[[1]]$tip.label)
write.csv(tip_label, "output_tree_reconstructed_phylo_tip-label.csv", row.names=FALSE, quote=FALSE) 

seqs <- read.csv("after_fitch_tree_seqs.csv")
name_lst <- c()
for (i in 1:length(seqs$id)){
  name <- seqs$id[i]
  name_subs <- substring(name,1,4)
  if (name_subs == 'NODE'){
    name_lst <- append(name_lst, name)
  }
}

frst_1625 <- tip_label$phylo_obj..1...tip.label
name_lst <- append(frst_1625, name_lst)
index <- 1:2445

ids <- data.frame(matrix(ncol=0, nrow=2445))
ids$number <- index
ids$name <- name_lst

write.csv(ids, "output_reconstructed_name_nums.csv", row.names=FALSE, quote=FALSE) 

