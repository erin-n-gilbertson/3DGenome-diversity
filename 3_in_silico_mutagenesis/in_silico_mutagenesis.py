# Colin M. Brand, University of California San Francisco, 09/19/2022
# modified from script by Evonne McArthur

import argparse
import itertools
import json
import numpy as np
import os
import pandas as pd
import pysam

from basenji import dataset, dna_io, seqnn
from scipy import stats

def parse_args():
	parser = argparse.ArgumentParser()
	
	parser.add_argument(
		"--fasta", type=str, required=True,
		help="Path to input FASTA.")

	parser.add_argument(
		"--input", type=str, required=True,
		help="Path to input file with variants. Formatted as a tab-delimited file with" 
		"chromosome, position, reference allele, alternate allele, and window start.")
		
	parser.add_argument(
		"--start", type=int, required=False, default=0,
		help="0-based line in input file at which to start in silico mutagenesis.")
		
	parser.add_argument(
		"--end", type=int, required=False,
		help="0-based line in input file at which to end in silico mutagenesis.")
		
	parser.add_argument(
		"--out", type=str, required=True, help="Path to output file.")
		
	args = parser.parse_args()
	return args

def main():
	args = parse_args()
	
	loadAkita()
	
	df = []
	
	with open(f'{args.input}', 'r') as input, open(f'{args.out}', 'w') as out:
		if args.end is None:
			file_end = len(input.readlines())
			print(file_end)
		elif args.end is not None:
			file_end = args.end
			
		for line in itertools.islice(input, args.start, file_end):
			line = line.split('\t')
			
			global chr
			chr = line[0]
			
			global variant_pos
			variant_pos = line[1]
			
			ref_allele = line[2]
			
			global alt_allele
			alt_allele = line[3]
			
			global window_start
			window_start = line[4]
			
			reference_prediction, alternate_prediction = get_contact_maps()
			mse, pearson, spearman = compare3Dmaps(reference_prediction, alternate_prediction)
			df.append({'chr':chr,'pos':variant_pos,'window':int(window_start),'ref':ref_allele,'alt':alt_allele,'mse':mse,'1-pearson':pearson,'1-spearman':spearman})
			print(f'{chr}: {variant_pos} completed')
	
	df = pd.DataFrame(df)[['chr','pos','window','ref','alt','mse','1-pearson','1-spearman']]
	df.to_csv(args.out, sep='\t', index=False)

def loadAkita():
	os.environ["CUDA_VISIBLE_DEVICES"] = '-1'

	import tensorflow as tf
	if tf.__version__[0] == '1':
		tf.compat.v1.enable_eager_execution()
		
	with open('../model/params.json') as params_file:
		params = json.load(params_file)
		params_model = params['model']
		params_train = params['train']

	global seqnn_model
	seqnn_model = seqnn.SeqNN(params_model)
	seqnn_model.restore('../model/model_best.h5')
	
def runAkitaPreds(seq):
	if len(seq) != 2**20: raise ValueError('len(seq) != seq_length')
	seq_1hot = dna_io.dna_1hot(seq)
	test_pred_from_seq = seqnn_model.model.predict(np.expand_dims(seq_1hot,0))
	return test_pred_from_seq

def get_contact_maps():
	args = parse_args()
	ref_fasta_open = pysam.Fastafile(f'{args.fasta}')
	ref_seq = ref_fasta_open.fetch(chr, int(window_start), int(window_start)+2**20).upper()
	ref_pred = runAkitaPreds(ref_seq)
	ref_pred_hff = ref_pred[:,:,0][0]

	pos_in_window = (int(variant_pos)-1)-int(window_start)
	alt_seq = ref_seq[0:pos_in_window] + alt_allele + ref_seq[pos_in_window+1:]
	alt_pred = runAkitaPreds(alt_seq)
	alt_pred_hff = alt_pred[:,:,0][0]
	return ref_pred_hff, alt_pred_hff

def compare3Dmaps(refpred, altpred):
    mse = np.mean(np.square(refpred  - altpred ))
    pearson = 1 - stats.pearsonr(refpred ,altpred )[0]
    spearman = 1 - stats.spearmanr(refpred ,altpred )[0]
    return mse, pearson, spearman

if __name__ == '__main__':
    main()
