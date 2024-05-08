#!/usr/bin/env python

# Will calculate pairwise sequence similarities across the genome in sliding 1 Mb windows for pairs of individual
# Input: takes two individual ids  (example: python runSeqComparisons.py "vindija AFR_ESN_female_HG03105")
# Output: returns a file with sequence comparisons and coverage between the individuals

print("Started python",flush=True)
import numpy as np
import pysam
import argparse
import pandas as pd
import os

BASE_PATH = '/wynton/group/capra/projects/modern_human_3Dgenome' # base directory level

#Wynton
BIN_PATH = os.path.join(BASE_PATH, "bin")  # where my scripts live
DATA_PATH = os.path.join(BASE_PATH, "data")  # where I dump new data 
RESULTS_PATH = os.path.join(BASE_PATH, "results")  # where I analyze results

SRC_PATH = os.path.join(BASE_PATH, "src")  # where any packages needed to run analyses live. I haven't started structuring things this way yet. 

COMP_PATH = os.path.join(DATA_PATH,"pairwise/hsmrca")
COMP_PATH = os.path.join(DATA_PATH,"pairwise/reference")

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--indivs", type=str, required=True,
        help="1KG individual to use")

    
    parser.add_argument(
        "--window_size_exp", type=str, required=True,
        help="individual 2")
    args = parser.parse_args()
    return args

### find file location for the individuals considered ###

def find_fastaFiles(indiv, chrm):
    fasta_dir = '%s/genomes' % DATA_PATH
    if (indiv == 'hg38_reference') | (indiv == 'GAGP_ancestral'):
        in_file_loc = '%s/baselines/%s.fa' % (fasta_dir,indiv)
    
    elif indiv == 'hsmrca_ancestral':
        in_file_loc = '%s/human_archaic_ancestor/human_archaic_ancestor_in_hg38_%s.fasta' % (fasta_dir,chrm)

    else:
        pop = indiv.split('_')[0]
        idi = indiv.split('_')[3]
        in_file_loc = '%s/1KG/%s/%s/%s_%s_hg38_full.fa' % (fasta_dir,pop, indiv,chrm,idi)

    indiv_fasta_open = pysam.Fastafile(in_file_loc)
    print("find fasta")
    print(indiv)
    print(in_file_loc)
    return indiv_fasta_open


        
#####################################






def main():
    args = parse_args()
    # input individuals
    indivname1=args.indivs[0]
    indivname2=args.indivs[1]


    print("Indiv1 = %s, Indiv2 = %s" % (indivname1, indivname2),flush=True)

    out_path = '%s/pairwise/sequence/scaling' % DATA_PATH

    windows = pd.read_table('%s/window_scale/windows%s.txt' % (DATA_PATH, args.window_size_exp))
    seqcomps = []
    for w in windows.index:
        chrm = windows.loc[w]['chr']
        start_loc = windows.loc[w]['sub_start']
        end_loc = windows.loc[w]['sub_end']
        try:
            indiv1_fasta_open = find_fastaFiles(indivname1, chrm)
            indiv2_fasta_open = find_fastaFiles(indivname2, chrm)

            print("starting sequence comparisons on %s" % start_loc,flush=True)
            indiv1_seq = indiv1_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()
            indiv2_seq = indiv2_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()
            print("length of seq1: %s, length of seq2: %s" % (str(len(indiv1_seq)), str(len(indiv2_seq))))
            # calculate coverage
            indiv1_coverage = np.mean([ 0 if s == "N" else 1 for s in indiv1_seq])
            indiv2_coverage = np.mean([ 0 if s == "N" else 1 for s in indiv2_seq])

            if len(indiv1_seq) !=0 :
                seqComp_raw = sum([1 if i1 == i2 else 0 for i1,i2 in zip(indiv1_seq,indiv2_seq)])/len(indiv1_seq)
            else:
                seqComp_raw = 'na'
            seqcomps.append(seqComp_raw)
        except OSError:
            print("Failed on chr: %s" % chrm)
            continue

    windows['seqComp'] = seqcomps
    windows.to_csv("%s/SeqComps%s_%s_vs_%s.txt" % (args.window_size_exp, out_path,indivname1,indivname2))

    return

if __name__ == '__main__':
    main()