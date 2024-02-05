# Erin N. Gilbertson, University of California San Francisco, 02/01/2024
# built on code from Evonne McArthur, David Rinker, and Colin Brand

import argparse
import numpy as np
from scipy import stats
from basenji import dataset, dna_io, seqnn
from scipy import stats
import pandas as pd

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
    indivs = pd.read_csv('%s/reference/1KG_unrelated_indivs.txt' % DATA_PATH, index_col=0)
    outfile = open(f'{args.chromosome}_{args.window}_comparisons.txt', 'w')

    for i in indivs['1KG']:
        i_fasta = pysam.Fastafile('/wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/1KG/%s/%s/%s_%s_hg38_full.fa' % (i.split('_')[0], i, args.chromosome, i.split('_')[-1]))
        i_seq = i_fasta.fetch(args.chromosome, args.window, args.window+2**20).upper()
        i_pred = runAkitaPreds(i_seq)
        i_pred = i_pred[:,:,0][0]
        i_fasta.close()
        
        for j in indivs['1KG']:
            if i != j:
                j_fasta = pysam.Fastafile('/wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/1KG/%s/%s/%s_%s_hg38_full.fa' % (j.split('_')[0], j, args.chromosome, j.split('_')[-1]))
                j_seq = j_fasta.fetch(args.chromosome, args.window, args.window+2**20).upper()
                j_fasta.close()
                j_pred = runAkitaPreds(j_seq)
                j_pred = j_pred[:,:,0][0]
                mse, spearman = comparePreds(i_pred, j_pred)
                outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(i, j, args.chromosome, args.window, mse, spearman))

    outfile.close()

if __name__ == '__main__':
	main()