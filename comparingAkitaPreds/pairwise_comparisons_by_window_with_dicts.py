# Erin N. Gilbertson, University of California San Francisco, 02/01/2024
# built on code from Evonne McArthur, David Rinker, and Colin Brand

import argparse
import numpy as np
from scipy import stats
from basenji import dataset, dna_io, seqnn
from scipy import stats
import pandas as pd
import os
from os import path
import json
import pysam

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('--chromosome', type = str, required = True, help = 'chromosome of window to compare')
	
	parser.add_argument('--window', type = int, required = True, help = 'start position of window to compare')
		
	args = parser.parse_args()
	return args


def loadAkita():
    os.environ["CUDA_VISIBLE_DEVICES"] = '-1'

    import tensorflow as tf
    if tf.__version__[0] == '1':
        tf.compat.v1.enable_eager_execution()

    with open('/wynton/group/capra/projects/modern_human_3Dgenome/bin/basenji/manuscripts/akita/params.json') as params_file:
        params = json.load(params_file)
        params_model = params['model']
        params_train = params['train']
        
    global seqnn_model
    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore('/wynton/group/capra/projects/modern_human_3Dgenome/bin/basenji/manuscripts/akita/model_best.h5')

def runAkitaPreds(seq):
    if len(seq) != 2**20: raise ValueError('len(seq) != seq_length')
    seq_1hot = dna_io.dna_1hot(seq)
    test_pred_from_seq = seqnn_model.model.predict(np.expand_dims(seq_1hot,0))
    return test_pred_from_seq

def comparePreds(individual_1, individual_2):
    mse = np.mean(np.square(individual_1 - individual_2)) 
    spearman = stats.spearmanr(individual_1, individual_2)[0]
    return (mse, spearman)


def main():
    loadAkita()
    args = parse_args()
    indivs = pd.read_csv('/wynton/group/capra/projects/modern_human_3Dgenome/data/reference/1KG_unrelated_indivs.txt', index_col=0)
    preds = {}
    for i in indivs['1KG'][:10]:
        if path.exists('/wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/1KG/%s/%s/%s_%s_hg38_full.fa' % (i.split('_')[0], i, args.chromosome, i.split('_')[-1])):
            i_fasta = pysam.Fastafile('/wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/1KG/%s/%s/%s_%s_hg38_full.fa' % (i.split('_')[0], i, args.chromosome, i.split('_')[-1]))
            i_seq = i_fasta.fetch(args.chromosome, args.window, args.window+2**20).upper()
            i_pred = runAkitaPreds(i_seq)
            i_pred = i_pred[:,:,0][0]
            i_fasta.close()
            preds[i] = i_pred

    comps = {}
    for i in preds.keys():
        for j in preds.keys():
            if (i != j) and ((i,j) not in comps.keys()):
                print(i,j) 
                mse, spearman = comparePreds(preds[i], preds[j])
                comps[(i,j)] = spearman
    df = pd.Series(comps).rename_axis(['Col1', 'Col2']).reset_index(name='spearman')
    df['divergence'] = 1- df['spearman']
    df = df.drop(columns=['spearman'])
    print(df.shape)
    triu = df.pivot(index='Col1', columns='Col2', values='divergence')
    tril = df.pivot(index='Col2', columns='Col1', values='divergence')
    sym = triu.fillna(tril)
    print(sym.shape)
    sym.to_csv('/wynton/group/capra/projects/modern_human_3Dgenome/data/pairwise/divergent_windows/%s_%s_comparisons_symmetric.txt' % (args.chromosome, args.window), sep='\t', index=False)
    df.to_csv('/wynton/group/capra/projects/modern_human_3Dgenome/data/pairwise/divergent_windows/%s_%s_comparisons_long.txt' % (args.chromosome, args.window), sep='\t', index=False)


if __name__ == '__main__':
	main()