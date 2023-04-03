#!/usr/bin/env python

import os
import sys
import json
import scipy.stats as stats
import pysam
import json
import re
import random

import numpy as np
import pandas as pd
from basenji import dataset, dna_io, seqnn
import tensorflow as tf
import pickle

np.random.seed(1337)

BASE_PATH = '/wynton/group/capra/projects/modern_human_3Dgenome'

BIN_PATH = os.path.join(BASE_PATH, "bin")  # where my scripts live
DATA_PATH = os.path.join(BASE_PATH, "data")  # where I dump new data 
RESULTS_PATH = os.path.join(BASE_PATH, "results")  # where I analyze results

FASTA_PATH = os.path.join(DATA_PATH, "genomes")

OUT_DIR = os.path.join(RESULTS_PATH, "expectedDist")

os.environ["CUDA_VISIBLE_DEVICES"] = '-1' ### run on CPU
print(tf.__version__)
if tf.__version__[0] == '1':
    tf.compat.v1.enable_eager_execution()

MODEL_DIR = '%s/basenji/manuscripts/akita/' % BIN_PATH
AKITA_PARAMS = 'params.json'
MODEL = 'model_best.h5'
TARGETS = "targets.txt"
STATS = "statistics.json"
#Setting up Akita model
### load params, specify model ###
print('load params and model')
params_file = os.path.join(MODEL_DIR, AKITA_PARAMS)
model_file  = os.path.join(MODEL_DIR, MODEL)
with open(params_file) as params_open:
    params = json.load(params_open)
    params_model = params['model']
    params_train = params['train']

seqnn_model = seqnn.SeqNN(params_model)

### restore model ###
print('restore model')
seqnn_model.restore(model_file)

print('read data paramenters')
data_stats_file = os.path.join(MODEL_DIR, "data/%s" % STATS)
with open(data_stats_file) as data_stats_open:
    data_stats = json.load(data_stats_open)
seq_length = data_stats['seq_length']
target_length = data_stats['target_length']
hic_diags =  data_stats['diagonal_offset']
target_crop = data_stats['crop_bp'] // data_stats['pool_width']
target_length1 = data_stats['seq_length'] // data_stats['pool_width']


target_length1_cropped = target_length1 - 2*target_crop
print('flattened representation length:', target_length)
print('symmetrix matrix size:', '('+str(target_length1_cropped)+','+str(target_length1_cropped)+')')

def main():
    mod_indiv = sys.argv[1]
    chrm = sys.argv[2]
    oneChromOneInd(mod_indiv, chrm)

    return


def runAkitaPreds(seq):
    print('run predictions')
    if len(seq) != 2**20: raise ValueError('len(seq) != seq_length')
    seq_1hot = dna_io.dna_1hot(seq)
    test_pred_from_seq = seqnn_model.model.predict(np.expand_dims(seq_1hot,0))
    return test_pred_from_seq

def find_inFileLoc(indiv, chrm):
    if indiv == 'hsmrca_ancestral':
        seq_dir = os.path.join(FASTA_PATH, "human_archaic_ancestor")
        file_path = os.path.join(seq_dir, ("human_archaic_ancestor_in_hg38_%s.fasta" % chrm))
    else:
        seq_dir = os.path.join(FASTA_PATH, "1KG")
        pop = indiv.split('_')[0]
        iid = indiv.split('_')[-1]
        file_path = os.path.join(seq_dir, ("%s/%s/%s_%s_hg38_full.fa" % (pop, indiv, chrm, iid)))
    return file_path


def oneChromOneInd(mod_indiv, chrm):

    chunks = pickle.load( open( "%s/reference/genome_chunks_dict.p" % DATA_PATH, "rb" ) )
    anc_indiv = 'hsmrca_ancestral'
    in_file_loc_anc_indiv = find_inFileLoc(anc_indiv, chrm)
    in_file_loc_mod_indiv = find_inFileLoc(mod_indiv, chrm)

    anc_indiv_fasta_open = pysam.Fastafile(in_file_loc_anc_indiv)
    mod_indiv_fasta_open = pysam.Fastafile(in_file_loc_mod_indiv)

    f = open("%s/%s_%s_empiricDist_test.tsv" % (OUT_DIR, mod_indiv,chrm),"w")

    for start_loc in chunks[chrm]:
        anc_indiv_seq = anc_indiv_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()
        mod_indiv_seq = mod_indiv_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()
        # Calculate the number of differences between the ancestral and modern sequences
        diff = [1 if anc != mod else 0 for anc,mod in zip(anc_indiv_seq,mod_indiv_seq)]
        diff_num = sum(diff)
        anc_pred  = runAkitaPreds(anc_indiv_seq)[:,:,0][0]
        
        diffs = dict()

        for i in range(len(anc_indiv_seq)):
            if diff[i] == 1: # make sure a difference is present
                if i == 0: # If the base is at the start of the sequence, the trinucleotide context is defined as the base itself and the following base
                    if (anc_indiv_seq[i+1] == 'N'):
                        triNuc = anc_indiv_seq[i]
                    else:
                        triNuc = "".join([anc_indiv_seq[i],anc_indiv_seq[i+1]])
                elif i == len(anc_indiv_seq)-1: #  # If the base is at the end of the sequence, the trinucleotide context is defined as the base itself and the preceding base
                    if (anc_indiv_seq[i-1] == "N"):
                        triNuc = anc_indiv_seq[i]
                    else:
                        triNuc = "".join([anc_indiv_seq[i-1],anc_indiv_seq[i]])
                else: # Otherwise, the trinucleotide context is defined as the three nucleotides immediately upstream and downstream of the base
                    if (anc_indiv_seq[i-1] == 'N') and (anc_indiv_seq[i+1] == 'N'):
                        triNuc = anc_indiv_seq[i]
                    elif (anc_indiv_seq[i-1] == 'N'):
                        triNuc = "".join([anc_indiv_seq[i],anc_indiv_seq[i+1]])
                    elif (anc_indiv_seq[i+1] == 'N'):
                        triNuc = "".join([anc_indiv_seq[i-1],anc_indiv_seq[i]])
                    else:
                        triNuc = "".join([anc_indiv_seq[i-1],anc_indiv_seq[i],anc_indiv_seq[i+1]])
                        
                # Add the archaic sequence at this position to the dictionary for this trinucleotide context
                if triNuc not in diffs:
                    diffs[triNuc] = [mod_indiv_seq[i]]
                else:     # If the trinucleotide context is already in the dictionary, append the archaic sequence at this position to the list for this trinucleotide context
                    diffs[triNuc].append(mod_indiv_seq[i])

        empirical_dist = []
        # Repeat the following process 100 times to create a null distribution of differences
        for iter in range(100): 
            print(iter,flush=True)
            # Create a randomly mutated version of the African individual
            random_mutated = anc_indiv_seq
            anc_ignore = anc_indiv_seq

            for k,v in diffs.items(): # For each of the differences
                    indexes = [m.start() for m in re.finditer(k, anc_ignore)] #find the sequence
                    indexes = random.sample(indexes,len(v)) #randomly sample
                    for i in indexes:
                        anc_ignore = anc_ignore[:i+1] + "N" + anc_ignore[i+2:] # mask the original sequence with an N so it will not try to resample this part of the sequence
                    for i,m in zip(indexes,v):
                        random_mutated = random_mutated[:i+1] + m + random_mutated[i+2:] # mutate trinucleotide in new posision
            
            random_pred  = runAkitaPreds(random_mutated)[:,:,0][0] #run akita on the shuffled sequence
            empirical_dist.append(1-stats.spearmanr(anc_pred, random_pred)[0]) # calculate the spearman on the shuffled sequence


        f.write(chrm + "\t" + str(start_loc) + "\t" + str(diff_num) + "\t" + "\t".join([str(x) for x in empirical_dist]) + "\n")
        f.flush()
    f.close()
    return

if __name__ == '__main__':
    main()