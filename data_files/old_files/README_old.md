# final-project
Code for CS 4775 Final Project!

## File Structure
* **data_files**: contains all data used, including FASTA sequences, newick tree files, phylogenetic tree data
* **ancestral_sequence_reconstruction**: contains file that was used to reconstruct the sequences at the internal nodes of the phylogenetic tree and save all required data
* **slac**: contains the file that parses the data after ancestral tree sequence reconstruction and applies the SLAC algorithm to this data, obtaining the required dN/dS values for each codon


## Pipeline and How to Run the Code
1. **Tree Construction**: 
   * This takes the aligned FASTA sequences from the *HA_cds.fa* file and infers a phylogeny from these. These sequences were downloaded from https://github.com/hzi-bifo/PatchDetection/tree/master/Software/Testdata/HA.
   * We then used the Generalized Time Reversible Model from the FastTree software. 
   
   * The FastTree software was first downloaded from http://www.microbesonline.org/fasttree/
   * In the terminal that opened up, the following command was typed: FastTree -gtr -nt < HA_cds.fa > tree_file
   * This allowed us to save the inferred phylogeny as *tree_file.nwk*.

2. **Ancestral Sequence Reconstruction**:
   * This utilizes the *HA_cds.fa* (which contains the sequences of the leaves of the tree) and *tree_file.nwk* (which contains the inferred phylogenetic tree) to reconstruct the sequences at the intermediate nodes (i.e., for each ancestor in the tree). 
   * The ACCTRAN version of the Fitch algorithm is used (which is an algorithm that tackles the small parsimony problem).
   * The file *fitch - final.R* can be run in RStudio, and will output a number of files:
   * *reconstructed_sequences_phyDat_randomly-sampled.csv* --> a list of the sequences corresponding to each node and leaf of the tree
   * *output_tree_reconstructed_phylo_edge.csv* --> a list of all the edges in the phylogenetic tree with the reconstructed sequences
   * *output_tree_reconstructed_phylo_edge-length.csv* --> a list of the lengths of all the edges in the phylogenetic tree
   * *output_tree_reconstructed_phylo_node-label.csv* --> the label of each leaf in the order the tree was traversed
   * *output_tree_reconstructed_phylo_tip-label.csv* --> the label of each leaf
   * *output_tree_reconstructed_phylo_nnode.csv* --> the Nnode variable of a phylo object in R

   * To run the code here, set the correct working directory as per your device in the code (in the line setwd(""))
   * Run the code in RStudio or VSCode with the R extension
   * The outputted data will be saved to your device. The correctly outputted files corresponding to *HA_cds.fa* and *tree_file.nwk* are in the data_files folder of this repository.

**OR**

   * The Fitch Algorithm is used (which is an algorithm that tackles the small parsimony problem).
   * This part of the code can be run in Python and utilizes the infer_ancestral_sequences method in the TreeAnc class of the TreeTime package. 
   * The infer_ancestral_sequences function outputs a tree with new labels and sequences corresponding to the intermediate nodes.
   * Thus, the outputted files from the code include:
   * *after_fitch_tree* --> a FASTA file representing the sequences for all nodes, including the intermediate nodes
   * *phy_tree.xml* --> a PhyloXML file representing the tree with the newly labeled intermediate nodes

   * To run the code here, the python file can simply be run in your IDE of choice.
   * If running the code in Jupyter Notebook, the following command can be called from the terminal: python final-project/ancestral_sequence_reconstruction/fitch.py
   * The outputted data will be saved to your device. The correctly outputted files corresponding to *HA_cds.fa* and *tree_file.nwk* are in the data_files folder of this repository.


3. **SLAC**:
   * This takes the files from the ancestral sequence reconstruction and constructs the phylogenetic tree using a Tree data structure defined in the python code.
   * It then runs the SLAC algorithm as defined by Suzuki and Gojobori (1999) on the tree and outputs dN/dS values for each codon of the sequence (barring codons that were removed due to a value of zero for the average number of synonymous sites at that codon site) and the determination of whether each site had negative, positive, or neutral selection.

   * To run the code here, the python file can simply be run in your IDE of choice.
   * If running the code in Jupyter Notebook, the following command can be called from the terminal: python final-project/slac/slac.py
   * Ensure that the required files follow the file structure as in this repository as they have been accessed in the python code using *final-project/data_files/{file name}*


## Dependencies (Packages Utilized)
For the Ancestral Sequence Reconstruction in R:
* ape
* phangorn
* seqinr
* phylotools
* TreeTools
* seqLogo

**OR**

* phylo-treetime
* BioPython

For the SLAC implementation in Python:
* numpy
* math
* itertools

  
## Sources
The sequences of the HA1 protein of the H3N2 Influenza virus used for this project were obtained from Klingen, T.R., Loers, J., Stanelle-Bertram, S. et al. Structures and functions linked to genome-wide adaptation of human influenza A viruses. Sci Rep 9, 6267 (2019).

The phylogenetic tree contruction was done using the FastTree algorithm (http://www.microbesonline.org/fasttree/).

**potentially remove following line**
The following website was referred to when conducting the ancestral sequence reconstruction: https://cran.r-project.org/web/packages/phangorn/vignettes/Ancestral.html.

The SLAC algorithm was based on the paper: Y Suzuki, T Gojobori. A method for detecting positive selection at single amino acid sites. Molecular Biology and Evolution. Volume 16, Issue 10. 1999. Pages 1315–1328.
