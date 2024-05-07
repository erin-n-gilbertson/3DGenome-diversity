print('python imports')
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle

import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist

from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist, jaccard
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import to_tree, ClusterNode, dendrogram
from typing import Dict, Tuple, List, Union, Optional
import os

from itertools import combinations

print('start python')
#local
BASE_PATH = "/".join(os.getcwd().split("/")) # base directory level

# DATA_PATH = os.path.join(BASE_PATH, "../../../downloads/")  # where I dump new data
# COMP_PATH = os.path.join(DATA_PATH,"1KGvs1KG")


#Wynton
BIN_PATH = os.path.join(BASE_PATH, "bin")  # where my scripts live
DATA_PATH = os.path.join(BASE_PATH, "data")  # where I dump new data 
RESULTS_PATH = os.path.join(BASE_PATH, "results")  # where I analyze results

SRC_PATH = os.path.join(BASE_PATH, "src")  # where any packages needed to run analyses live. I haven't started structuring things this way yet. 

COMP_PATH = os.path.join(DATA_PATH, "pairwise/all_windows")
# COMP_PATH = os.path.join(DATA_PATH,"pairwise/hsmrca")
# COMP_PATH = os.path.join(DATA_PATH,"pairwise/reference")



# Create pickle data structures 3D by window
print('setup data structures')
windows = pd.read_table('%s/intermediates/windows_to_keep.csv' % DATA_PATH, sep=',', index_col=[1,2]).drop(columns=['Unnamed: 0'])
indivs = pd.read_csv('/wynton/group/capra/projects/modern_human_3Dgenome/data/reference/1KG_unrelated_indivs.txt', index_col=0)
comp_list = list(combinations(list(indivs['1KG']), 2))
df = pd.read_table('%s/chr9_9961472_comparisons_dict.txt' % COMP_PATH)

count = 0
dict_3d = {} 
matrix_3d = pd.DataFrame(df[['indiv1','indiv2']])
print('start looping')
for w in list(windows.index):
    print(w)
    try:
        df = pd.read_table('%s/%s_%s_comparisons_dict.txt' % (COMP_PATH, w[0], w[1]))
        df.columns = ['indiv1','indiv2',w]
        matrix_3d = pd.concat([matrix_3d, df[w]],axis=1)
        count +=1 

    except:
        print("error on %s" % w)


    print(count)
    triu = df.pivot(index='indiv1', columns='indiv2', values=w)
    tril = df.pivot(index='indiv2', columns='indiv1', values=w)
    sym = triu.fillna(tril)
    new_col = sym.loc['AFR_ACB_female_HG01880'].T
    new_col = pd.concat([pd.Series([np.nan], index=[new_col.name]), new_col])
    new_row = sym['SAS_STU_male_HG04229'].T
    new_row = pd.concat([new_row, pd.Series([np.nan], index=[new_row.name])])
    sym['AFR_ACB_female_HG01880'] = np.nan
    sym = pd.concat([sym['AFR_ACB_female_HG01880'], sym[sym.columns[:-1]]], axis=1)
    sym.loc['SAS_STU_male_HG04229'] = new_row.values
    sym['AFR_ACB_female_HG01880'] = new_col.values
    dict_3d[w] = sym

pickle.dump( dict_3d, open( "%s/dict_3d_all_pairs.p" % COMP_PATH, "wb" ) )
pickle.dump( matrix_3d, open( "%s/matrix_3d_all_pairs.p" % COMP_PATH, "wb" ) )


