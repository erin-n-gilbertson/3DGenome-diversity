print('python imports')
import pandas as pd
import numpy as np

import pickle

import os
import argparse

print('start python')

BASE_PATH = '/wynton/group/capra/projects/modern_human_3Dgenome' # base directory level
#Wynton
BIN_PATH = os.path.join(BASE_PATH, "bin")  # where my scripts live
DATA_PATH = os.path.join(BASE_PATH, "data")  # where I dump new data 
RESULTS_PATH = os.path.join(BASE_PATH, "results")  # where I analyze results

COMP_PATH = os.path.join(DATA_PATH, "pairwise/all_windows")

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--chromosome", type=str, required=True,
        help="chromosome")

    args = parser.parse_args()
    return args










def main():
    args = parse_args()
    # Create pickle data structures 3D by window
    print('setup data structures')
    windows = pd.read_table('%s/intermediates/windows_to_keep.csv' % DATA_PATH, sep=',', index_col=0)
    #indivs = pd.read_csv('/wynton/group/capra/projects/modern_human_3Dgenome/data/reference/1KG_unrelated_indivs.txt', index_col=0)
    #comp_list = list(combinations(list(indivs['1KG']), 2))
    df = pd.read_table('%s/chr9_9961472_comparisons_dict.txt' % COMP_PATH)

    windows = windows[windows.chr==args.chromosome]

    count = 0
    dict_3d = {} 
    matrix_3d = pd.DataFrame(df[['indiv1','indiv2']])

    print('start looping')
    for w in list(windows.index):
        print(w)
        wchr = windows.loc[w]['chr']
        wpos = windows.loc[w]['windowStartPos']
        try:
            df = pd.read_table('%s/%s_%s_comparisons_dict.txt' % (COMP_PATH, wchr, wpos))
            df.columns = ['indiv1','indiv2',(wchr, wpos)]
            matrix_3d = pd.concat([matrix_3d, df[w]],axis=1)
            count +=1 

        except:
            print("error on %s, %s" % (wchr, wpos))


        print(count)
        triu = df.pivot(index='indiv1', columns='indiv2', values=[(wchr, wpos)])
        tril = df.pivot(index='indiv2', columns='indiv1', values=[(wchr, wpos)])
        sym = triu.fillna(tril)
        new_col = sym.loc['AFR_ACB_female_HG01880'].T
        new_col = pd.concat([pd.Series([np.nan], index=[new_col.name]), new_col])
        new_row = sym['SAS_STU_male_HG04229'].T
        new_row = pd.concat([new_row, pd.Series([np.nan], index=[new_row.name])])
        sym['AFR_ACB_female_HG01880'] = np.nan
        sym = pd.concat([sym['AFR_ACB_female_HG01880'], sym[sym.columns[:-1]]], axis=1)
        sym.loc['SAS_STU_male_HG04229'] = new_row.values
        sym['AFR_ACB_female_HG01880'] = new_col.values
        dict_3d[(wchr, wpos)] = sym

    pickle.dump( dict_3d, open( "%s/dict_3d_all_pairs_%s.p" % (COMP_PATH, args.chromosome), "wb" ) )
    pickle.dump( matrix_3d, open( "%s/matrix_3d_all_pairs_%s.p" % (COMP_PATH, args.chromosome), "wb" ) )  
    return

if __name__ == '__main__':
    main()
