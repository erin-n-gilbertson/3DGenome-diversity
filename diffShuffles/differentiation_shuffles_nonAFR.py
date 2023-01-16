import pandas as pd
import os

import pickle
import random
import numpy as np 
import pybedtools

def calc_one_pop_stats(within_comps, between_comps):
    
    tmp = pd.DataFrame(index=['within','between','all'], columns = comp_list_div.columns[2:])
    tmp.loc['within'] = comp_list_div.loc[within_comps].mean()
    tmp.loc['between'] = comp_list_div.loc[between_comps].mean()
    tmp.loc['all'] = comp_list_div.loc[all_comps].mean()
    tmp.loc['new_fst'] = (tmp.loc['all'] - tmp.loc['within']) / tmp.loc['all']
    return tmp.T


comp_list=pd.read_csv('/wynton/group/capra/projects/modern_human_3Dgenome/results/comp_tables/pairwise_subsample_genomewide_averages.csv', index_col=0)
# Only keep when excluding AFR
comp_list = comp_list[(comp_list.super1 != 'AFR') & (comp_list.super2 != 'AFR')]

comp_list_div = pd.read_csv('/wynton/group/capra/projects/modern_human_3Dgenome/results/comp_tables/pairwise_subsample_divergence_per_region.csv',index_col=0, header=[0,1])


all_comps = list(comp_list.index)

num_shuffles = 10000

shuffle_fsts = pd.DataFrame(index=comp_list_div.columns[2:], columns = range(1,num_shuffles+1))
#shuffle_diffs = s
indivs = list(np.unique(np.concatenate([comp_list.ind1,comp_list.ind2])))

for i in range(1,num_shuffles+1):
    if (i < 100) or (i%100==0):
        print(i)
    shuffle = random.sample(indivs, len(indivs))
    within_inds = shuffle[:26]
    within_comps = list(comp_list[(comp_list.ind1.isin(within_inds)) & (comp_list.ind2.isin(within_inds))].index)
    between_comps = list(comp_list.index[comp_list.index.isin(within_comps)==False])
    
    diff = calc_one_pop_stats(within_comps, between_comps)
    shuffle_fsts[i] = diff['new_fst']
    if i%500==0 :
        shuffle_fsts.to_csv('/wynton/group/capra/projects/modern_human_3Dgenome/results/differentiation_shuffles/basic_%s_shuffles_differentiation_nonAFR.csv' % i)


