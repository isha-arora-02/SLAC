from treetime import TreeAnc
from Bio import Phylo
from Bio import Align
from Bio import AlignIO
from Bio.Phylo import PhyloXMLIO
import csv

tree = Phylo.read('final-project-main/data_files/tree_file.nwk', 'newick')
phy_seq = AlignIO.read('final-project-main/data_files/HA_cds.fa', 'fasta')

# create a TreeAnc object
phy_tree = TreeAnc(tree=tree, aln=phy_seq)

# infer ancestral sequences on tree
phy_tree.infer_ancestral_sequences('fitch')

# obtain sequences of all nodes including ancestral sequences
after_fitch_tree = phy_tree.get_tree_dict()

# write all sequences and ancestral sequences to FASTA file
AlignIO.write(after_fitch_tree, 'final-project-main/data_files/after_fitch_tree_seqs.fasta', 'fasta')

# obtain phylogenetic tree as an PhyXML object
Phylo.write(phy_tree.tree, 'final-project-main/data_files/phy_tree.xml', 'phyloxml')
