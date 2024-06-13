print('python imports')

import pandas as pd
import numpy as np
import pickle 
import scipy.cluster.hierarchy as sch
import scipy.cluster.hierarchy as hc
from scipy.cluster.hierarchy import to_tree, ClusterNode, dendrogram
from typing import Dict, Tuple, List, Union, Optional
from ete3 import Tree
import os


print('start python')

BASE_PATH = '/wynton/group/capra/projects/modern_human_3Dgenome' # base directory level
#Wynton
BIN_PATH = os.path.join(BASE_PATH, "bin")  # where my scripts live
DATA_PATH = os.path.join(BASE_PATH, "data")  # where I dump new data 
RESULTS_PATH = os.path.join(BASE_PATH, "results")  # where I analyze results

COMP_PATH = os.path.join(DATA_PATH, "pairwise/all_windows")


def make_dict_3d():
    print('make dict 3d')
    dict_3d = {}
    slices = [(1,5), (5,10), (10,15), (15,20), (20,23)]

    for i in slices:
        print(i)
        dicti = pickle.load( open( "%s/dict_3d_all_pairs_%s-%s.p" % (COMP_PATH, i[0], i[1]), "rb" ) )
        dict_3d.update(dicti)

    print('save dict 3d')

    pickle.dump( dict_3d, open( "%s/dict_3d_all_pairs.p" % (COMP_PATH), "wb" ) )

    print('done make dict 3d')

    return dict_3d

def one_window_tree(dict3d, w, tree_summary, tree_complete, link_method, idx):
        window_df = dict_3d[w].loc[idx][idx]
        length = len(idx)
        array = window_df.reindex(index=idx, columns=idx).fillna(0, downcast='infer').to_numpy()
        condensed = array[np.triu_indices(length, k = 1)]
        Z = sch.linkage(condensed, method = link_method)
        dendrogram = sch.dendrogram(Z, labels=idx)
        whole_tree = [Z, dendrogram]

        cluster_sample_IDs = dendrogram['leaves']
        cluster_IDs = dendrogram['leaves_color_list']
        clusters_dict = dict(zip(cluster_sample_IDs, cluster_IDs))

        cluster_IDs_list = []
        for key, value in sorted(clusters_dict.items()):
                cluster_IDs_list.append(value)

        top_tree_y = dendrogram['dcoord'][-1]
        window_stats = cluster_IDs_list + top_tree_y


        return window_stats, whole_tree

def _scipy_tree_to_newick_list(node: ClusterNode, newick: List[str], parentdist: float, leaf_names: List[str]) -> List[str]:
    """Construct Newick tree from SciPy hierarchical clustering ClusterNode

    This is a recursive function to help build a Newick output string from a scipy.cluster.hierarchy.to_tree input with
    user specified leaf node names.

    Notes:
        This function is meant to be used with `to_newick`

    Args:
        node (scipy.cluster.hierarchy.ClusterNode): Root node is output of scipy.cluster.hierarchy.to_tree from hierarchical clustering linkage matrix
        parentdist (float): Distance of parent node of `node`
        newick (list of string): Newick string output accumulator list which needs to be reversed and concatenated (i.e. `''.join(newick)`) for final output
        leaf_names (list of string): Leaf node names

    Returns:
        (list of string): Returns `newick` list of Newick output strings
    """
    if node.is_leaf():
        return newick + [f'{leaf_names[node.id]}:{parentdist - node.dist}']

    if len(newick) > 0:
        newick.append(f'):{parentdist - node.dist}')
    else:
        newick.append(');')
    newick = _scipy_tree_to_newick_list(node.get_left(), newick, node.dist, leaf_names)
    newick.append(',')
    newick = _scipy_tree_to_newick_list(node.get_right(), newick, node.dist, leaf_names)
    newick.append('(')
    return newick


def to_newick(tree: ClusterNode, leaf_names: List[str]) -> str:
    """Newick tree output string from SciPy hierarchical clustering tree

    Convert a SciPy ClusterNode tree to a Newick format string.
    Use scipy.cluster.hierarchy.to_tree on a hierarchical clustering linkage matrix to create the root ClusterNode for the `tree` input of this function.

    Args:
        tree (scipy.cluster.hierarchy.ClusterNode): Output of scipy.cluster.hierarchy.to_tree from hierarchical clustering linkage matrix
        leaf_names (list of string): Leaf node names

    Returns:
        (string): Newick output string
    """
    newick_list = _scipy_tree_to_newick_list(tree, [], tree.dist, leaf_names)
    return ''.join(newick_list[::-1])

def write_trees_to_file(trees, output_file):
    with open(output_file, "w") as f:
        for name, tree in trees.items():
            #f.write(">" + name + "\n") # Write tree name as a header
            f.write(tree.write(format=1).strip() + "\n") # Write tree in Newick format

def main():
    indivs = pd.read_csv('/wynton/group/capra/projects/modern_human_3Dgenome/data/reference/1KG_unrelated_indivs.txt', index_col=0)
    windows = pd.read_table('%s/intermediates/windows_to_keep.csv' % DATA_PATH, sep=',', index_col=[1,2]).drop(columns=['Unnamed: 0'])
    dict_3d = make_dict_3d()

    idx = list(indivs['1KG'])
    idx.remove('SAS_ITU_male_HG04060')
    sub_idx = []
    for i in idx:
        sub_idx.append(i[4:7])

    non_adm = []
    for i in range(len(sub_idx)):
        p = sub_idx[i]
        if p not in ['ACB','ASW','CLM','MXL','PEL','PUR']:
            non_adm.append(i)

    non_adm_ids = [idx[i] for i in non_adm]

    print('cluster into trees')
    # Non admixed version
    tree_summary = {}
    tree_complete = {}
    #takes about 3 minutes 20 seconds +/-
    for i in range(len(list(windows.index))):
        w = list(windows.index)[i]
        window_stats, whole_tree = one_window_tree(dict_3d, w, tree_summary, tree_complete, link_method='complete', idx=non_adm_ids)
        tree_summary[w] = window_stats
        tree_complete[w] = whole_tree

        if i%100==0:
            print(i)

    pickle.dump( tree_summary, open( "%s/tree_summary_nonadm.p" % (COMP_PATH), "wb" ) )
    pickle.dump( tree_complete, open( "%s/tree_complete_nonadm.p" % (COMP_PATH), "wb" ) )

    trees = {}
    for w in list(windows.index):
        Z = tree_complete[w][0]
        T = Tree(to_newick(hc.to_tree(Z), non_adm_ids))
        trees[w] = T

    output_file = "%s/trees_nonadm.txt" % COMP_PATH
    write_trees_to_file(trees, output_file)
    
    return


if __name__ == '__main__':
    main()
