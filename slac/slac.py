"""
Reimplementation of SLAC algorithm.
Authors:
Ori Ben Yossef [oby3]
Isha Arora [ia93]
"""

import argparse
import numpy as np
import math
from itertools import permutations
from statsmodels.stats.proportion import binom_test




class Tree:
    """
    Initialize a node of a tree.

    Arguments:
    name:                  the name of that node
    sequence:              the sequence at that node itself
    branch_lst:            a list of dictionaries consisting of the children of the current node and the length of the branch from
                           the intermediate node to each child and this list takes the format:
                           [{'child':____, 'length':___},...]
    parent:                backpointer to the parent node
                           (this argument is None if there is no parent node)
    """
    def __init__(self, name, sequence, branch_lst, parent):
        self.name = name
        self.seq = sequence
        self.branches = branch_lst
        self.parent = parent



"""Dictionary of Codons and their corresponding Amino Acids"""
amino = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"$", "TAG":"$",
    "TGT":"C", "TGC":"C", "TGA":"$", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}


"""
Helper Function Used for Testing Parse_Tree().

Arguments:
tree:    an object of class Tree

Returns:
tree.name:                 the name of the root of tree
tree.seq:                  the sequence at that node itself
tree.branches:             a list of dictionaries consisting of the children of the current node and the length of the branch from
                           the intermediate node to each child
tree.parent:               backpointer to the parent node
"""
def data_of_tree(tree):
    return [tree.name, tree.seq[0:10], tree.branches, tree.parent]


"""
Parse the raw text file representing a tree,
and convert that to a data structure representing a tree.

Arguments: none

Returns:
tree:    an object of class Tree
"""
def parse_tree():

    # edge lengths
    with open("final-project-main/data_files/output_tree_reconstructed_phylo_edge-length.csv") as l_file:
        l_stream = l_file.read()
        l_all_lines = l_stream.split("\n")
        l_lines = l_all_lines[1:-1] # omit title and final empty string

    # edges | (parent number, child number) pairs
    with open("final-project-main/data_files/output_tree_reconstructed_phylo_edge.csv") as e_file:
        e_stream = e_file.read()
        e_all_lines = e_stream.split("\n")
        e_lines = e_all_lines[1:-1] # omit title and final empty string
        if len(l_lines) != len(e_lines):
            raise "different number of edges and edge lengths"
        e_lines_pairs = list(map( ( lambda x : x.split(",") ), e_lines ))
        n_edges = len(e_lines_pairs)

    # sequences | (name, sequence) pairs
    with open("final-project-main/data_files/after_fitch_tree_seqs.csv") as s_file:
        s_stream = s_file.read()
        s_all_lines = s_stream.split("\n")
        s_lines = s_all_lines[1:-1] # omit title and final empty string
        if len(s_lines) != len(e_lines) + 1:
            raise "incompatible number of edges and vertices"
        s_lines_pairs = list(map( ( lambda x : x.split(",") ), s_lines ))
        n_nodes_and_tips = len(s_lines_pairs)

    # names | (number, name) pairs
    with open("final-project-main/data_files/output_reconstructed_name_nums.csv") as n_file:
        n_stream = n_file.read()
        n_all_lines = n_stream.split("\n")
        n_lines = n_all_lines[1:-1] # omit title and final empty string
        if len(s_lines) != len(n_lines):
            raise "incompatible number of names and sequences"
        n_lines_pairs = list(map( ( lambda x : x.split(",") ), n_lines ))

    # root
    with open("final-project-main/data_files/after_fitch_root.csv") as r_file:
        r_stream = r_file.read()
        r_all_lines = r_stream.split("\n")
        r_first_line = r_all_lines[0] # omit final empty string
        r_first_line_cols = r_first_line.split(",") # separate 1st col from rest
        root_name = r_first_line_cols[0]

    # associate node id numbers with node names
    node_no = {}
    for name_num in n_lines_pairs:
        node_no[name_num[0]] = name_num[1]

    # initialize dict of trees. keyword is node name. value is Tree object.
    dict_of_trees = {}

    # link each child to branch length

    # first, construct all the trees.
    # recall constructor format: Tree(name, sequence, branch_lst, parent)
    # give each tree a sequence, but no parent or children.
    for name_seq in s_lines_pairs:
        dict_of_trees[name_seq[0]] = Tree(name_seq[0], name_seq[1], [], None)

    # then, go over every edge and assign the appropriate child/parent pointers
    for i_edge in range(n_edges):
        pnum_cnum = e_lines_pairs[i_edge]
        parent_name = node_no[pnum_cnum[0]]
        child_name = node_no[pnum_cnum[1]]
        child_distance = float(l_lines[i_edge])

        # assign child to parent
        dict_of_trees[parent_name].branches.append({'child': dict_of_trees[child_name], 'length': child_distance})
        #{'child':____, 'length':___}

        # assign parent to child
        dict_of_trees[child_name].parent = dict_of_trees[parent_name]

    return dict_of_trees[root_name]


"""
Helper function to calculate the expected synonymous and nonsynonymous changes given one codon.

Arguments:
codon:              a codon

Returns:
exp_syn:            number of expected synonymous changes
exp_nonsyn:         number of expected non-synonymous changes
"""
def exp_vals(codon):
    acgt = ['A', 'C', 'G', 'T']

    # make list of all possible codons with one mutation
    list_from_codon = list(codon)
    new_codons = []
    for ipos in range(3):
        for new_letter in acgt:
            list_from_new_codon = list_from_codon.copy()
            list_from_new_codon[ipos] = new_letter
            new_codon = ''.join(list_from_new_codon)
            if new_codon != codon:
                new_codons.append(new_codon)

    # count possible codons with expected synonymous, nonsynonymous mutation
    # (this will be divided by 3 in the count)
    exp_syn = 0
    exp_nonsyn = 0
    for new_codon in new_codons:
        if amino[new_codon] == amino[codon]:
            exp_syn += 1
        else:
            exp_nonsyn += 1
    exp_syn /= 3
    exp_nonsyn /= 3

    return exp_syn, exp_nonsyn


"""
Helper function to calculate the actual synonymous and nonsynonymous changes between two codons.

Arguments:
codon1:         the first codon
codon2:         the second codon that will be compared with the first codon

Returns:
avg_exp_syn:        average number of synonymous changes between the two codons
avg_exp_nonsyn:     average number of actual non-synonymous changes between the two codons
"""
def avg_exp_vals(codon1, codon2):

    # if the codons are equal
    if codon1 == codon2:
        # the exp_vals for both codons are equal, and thus the avg_exp_syn and avg_exp_nonsyn are equal to the exp_vals of either codon
        avg_exp_syn, avg_exp_nonsyn = exp_vals(codon1)
    else:

        # calculate the number of non-synonymous changes
        # average over all paths
        mutated_positions = []
        for ipos in range(3):
            if codon1[ipos] != codon2[ipos]:
                mutated_positions.append(ipos)
        paths = permutations(mutated_positions)
        acc_exp_syn, acc_exp_nonsyn = 0, 0 # running total over all paths
        for path in list(paths):
            path_exp_syn, path_exp_nonsyn = exp_vals(codon1) # running total over current path
            current_codon = list(codon1)
            for mutated_pos in path:
                new_codon = current_codon.copy()
                new_codon[mutated_pos] = codon2[mutated_pos]
                # ES and EN for intermediate codon in mutation path:
                intmd_exp_syn, intmd_exp_nonsyn = exp_vals(''.join(new_codon))
                # add it to running total
                path_exp_syn += intmd_exp_syn
                path_exp_nonsyn += intmd_exp_nonsyn
                current_codon = new_codon
            acc_exp_syn += path_exp_syn / (len(mutated_positions) + 1)
            acc_exp_nonsyn += path_exp_nonsyn / (len(mutated_positions) + 1)
        # divide to find average
        avg_exp_syn = acc_exp_syn / math.factorial(len(mutated_positions))
        avg_exp_nonsyn = acc_exp_nonsyn / math.factorial(len(mutated_positions))

    return avg_exp_syn, avg_exp_nonsyn



"""
Helper function to calculate the actual synonymous and nonsynonymous changes between two codons.

Arguments:
codon1:             the first codon
codon2:             the second codon that will be compared with the first codon

Returns:
actual_syn:         number of actual synonymous changes
actual_nonsyn:      number of actual non-synonymous changes
"""
def num_syn_nonsyn(codon1, codon2):

    # initializing values
    actual_syn = 0
    actual_nonsyn = 0
    diff_bases = 0

    # storing the amino acid corresponding to each codon
    aa1 = amino[codon1]
    aa2 = amino[codon2]

    # if the codons are the same
    if codon1 == codon2:
        return actual_syn, actual_nonsyn

    # if codons are not equal
    else:
        for ibase in range(3):
            if codon1[ibase] != codon2[ibase]:
                diff_bases += 1

        # if the resultant amino acids from the codons are equal, the changes are all synonymous
        if aa1 == aa2:
            actual_syn = diff_bases

        # if the resultant amino acids from the codons are NOT equal, at least one change is NON-synonymous
        else:
            # calculate the number of non-synonymous changes
            # average over all paths
            mutated_positions = []
            for ipos in range(3):
                if codon1[ipos] != codon2[ipos]:
                    mutated_positions.append(ipos)
            paths = permutations(mutated_positions)
            acc_nonsyn = 0 # running total over all paths
            for path in list(paths):
                nonsyn = 0
                current_codon = list(codon1)
                for mutated_pos in path:
                    new_codon = current_codon.copy()
                    new_codon[mutated_pos] = codon2[mutated_pos]

                    if amino[''.join(current_codon)] != amino[''.join(new_codon)]:
                        nonsyn += 1
                    current_codon = new_codon
                acc_nonsyn += nonsyn
            actual_nonsyn = acc_nonsyn / math.factorial(len(mutated_positions))
            actual_syn = len(mutated_positions) - actual_nonsyn

    return actual_syn, actual_nonsyn



""""
Recursive helper function where data is collected. Recurses through the tree by starting at the root
to calculate required data points for a particular codon of the gene at a global scale.

Arguments:
tree:      pointer to the root of the tree
start:     index of the start of a codon in the gene

Returns:
lst:       list with the required values for that particular codon and is of the form: [sum_tbl, sum_tbl_syn, sum_tbl_non, syns, nons] where
                - sum_tbl = Denominator of EN and ES = total branch lengths
                - sum_tbl_syn = Numerator of ES = total sum of (branch length * expected syn) for one codon
                - sum_tbl_non = Numerator of EN = total sum of (branch length * expected non-syn) for one codon
                - syns = NS = actual synonymous for one codon
                - nons = NN = actual non-synonymous for one codon
"""
def collect_data(tree, start):

    # initialize list that will store values to be returned
    lst = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

    # if the node is a leaf
    if (tree.branches == []):
        return lst

    # store the codon for which the data needs to be collected from
    # the sequence of the current node
    curr_codon = tree.seq[start:(start+3)]

    # looping over all children of the current node
    for child_dict in tree.branches:
        child = child_dict['child']
        branch_length = child_dict['length']

        # store the codon of the child node
        child_codon = child.seq[start:(start+3)]

        # average expected values of synonymous and nonsynonymous mutations between the current codon and the child codon
        avg_exp_syn_child, avg_exp_nonsyn_child = avg_exp_vals(curr_codon, child_codon)

        # update the values in the list to be returned as per lecture 25 notes, slides 29-30
        lst[0] += (branch_length)
        lst[1] += (avg_exp_syn_child * branch_length)
        lst[2] += (avg_exp_nonsyn_child * branch_length)

        # actual values of synonymous and nonsynonymous mutations between the intermediate codon and the child codon
        inter_child_syn, inter_child_nonsyn = num_syn_nonsyn(curr_codon, child_codon)

        # update the values in the list to be returned as per lecture 25 notes, slides 29-30
        lst[3] += (inter_child_syn)
        lst[4] += (inter_child_nonsyn)

        # recursive call on the child tree/node
        lst = lst + collect_data(child, start)

    return lst



"""
SLAC Algorithm.

Arguments:
tree:    an object of class Tree

Returns:
dnds:    dN/dS score for each codon in gene
p:       probability that a random mutation in the gene is synonymous
"""
def slac(tree):

    # if the tree is a single node
    if (tree.branches == []):
        return -1, -1

    # find root sequence
    curr_seq = tree.seq

    # Initialize empty np.arrays
    sum_tbl_lst = np.array([])
    sum_tbl_syn_lst = np.array([])
    sum_tbl_non_lst = np.array([])
    syns_lst = np.array([])
    nons_lst = np.array([])

    # traverse the sequence of the gene, determine the indices of codon, and calculate required values for each codon of the gene
    for icodon in range(len(tree.seq) // 3):
        start = icodon*3
        lst = collect_data(tree, start)

        sum_tbl = lst[0]
        sum_tbl_syn = lst[1]
        sum_tbl_non = lst[2]
        syns = lst[3]
        nons = lst[4]

        # append calculated values to respective lists
        sum_tbl_lst = np.append(sum_tbl_lst, sum_tbl)
        sum_tbl_syn_lst = np.append(sum_tbl_syn_lst, sum_tbl_syn)
        sum_tbl_non_lst = np.append(sum_tbl_non_lst, sum_tbl_non)
        syns_lst = np.append(syns_lst, syns)
        nons_lst = np.append(nons_lst, nons)

    # Calculate Expected Synonymous and Non-synonymous values for each Codon
    es = np.divide(sum_tbl_syn_lst, sum_tbl_lst)
    en = np.divide(sum_tbl_non_lst, sum_tbl_lst)

    # Remove codons with average value of zero in syns_lst and nons_lst from all data
    bool_lst_syns = [s != 0 for s in syns_lst]
    bool_lst_nons = [n != 0 for n in nons_lst]
    bool_lst = [(s and n) for s, n in zip(bool_lst_syns, bool_lst_nons)]

    lst_removed_codons = ['kept' if x == True else 'removed' for x in bool_lst]

    syns_lst = syns_lst[bool_lst]
    nons_lst = nons_lst[bool_lst]
    es = es[bool_lst]
    en = en[bool_lst]

    # Obtain list of indices of the codons that are remaining in the dataset after filtering
    l = len(bool_lst)
    indices = [x+1 if bool_lst[x] == True else None for x in list(range(l))]
    indices = [i for i in indices if i is not None]

    # Calculate dN, dS, and dN/dS for each Codon
    ds = np.divide(syns_lst, es)
    dn = np.divide(nons_lst, en)

    dnds = np.divide(dn, ds)

    # Conduct significance test as in paper
    sc = np.ndarray.round(syns_lst)
    nc = np.ndarray.round(nons_lst)
    expected_syn = np.divide(es, (es + en))
    expected_non = np.divide(en, (es + en))

    p_val_syn = []
    p_val_non = []
    for i in range(len(sc)):
        p_syn = binom_test(count=sc[i], nobs=(sc[i] + nc[i]), prop=expected_syn[i], alternative='larger')
        p_val_syn.append(p_syn)

        p_non = binom_test(count=nc[i], nobs=(sc[i] + nc[i]), prop=expected_non[i], alternative='larger')
        p_val_non.append(p_non)

    p_neg_sel = p_val_syn
    p_pos_sel = p_val_non

    dnds = list(zip(indices, dnds))
    p_neg_sel = list(zip(indices, p_neg_sel))
    p_pos_sel = list(zip(indices, p_pos_sel))

    # evaluate which codons have statistically significant selection
    c_neg_sel = []
    c_pos_sel = []
    for codon_pvalue in p_neg_sel:
        if codon_pvalue[1] < 0.05:
            c_neg_sel.append(codon_pvalue[0])
    for codon_pvalue in p_pos_sel:
        if codon_pvalue[1] < 0.05:
            c_pos_sel.append(codon_pvalue[0])


    return dnds, p_neg_sel, p_pos_sel, c_neg_sel, c_pos_sel


def main():

    tree = parse_tree()
    dnds, p_neg, p_pos, c_neg, c_pos = slac(tree)

    print('Indexing for codons starts at 1.')
    print('Alternative Hypothesis: there is negative selection at a codon. Requires p-value < 0.05 -->\n', p_neg)
    print('Alternative Hypothesis: there is positive selection at a codon. Requires p-value < 0.05 -->\n', p_pos)
    print('Negative selection was found in the following ' + str(len(c_neg)) + ' codon(s):\n', c_neg)
    print('Positive selection was found in the following ' + str(len(c_pos)) + ' codon(s):\n', c_pos)
    print('The dN/dS values for each codon of this gene are:\n', dnds)


if __name__ == "__main__":
    main()
