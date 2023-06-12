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
		"--fasta", type=str, required=True, help="Path to input FASTA.")

	parser.add_argument(
		"--input", type=str, required=True,
		help="Path to input file with variants. Formatted as a tab-delimited file with" 
		"chromosome, position, reference allele, alternate allele, window start, and individual.")
		
	args = parser.parse_args()
	return args

def main():
	args = parse_args()
	
	loadAkita()
	
	with open(f'{args.input}', 'r') as input:			
		for line in input:
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
			
			global individual
			individual = line[5].strip()
			
			alt_pred = alternate_map()
			
			with open(f'IDW_3d_modifying_variant_predictions/{chr}_{window_start}_{variant_pos}.txt', 'w') as prediction_out:
				prediction_out.write(chr + "\t" + str(window_start) + "\t" + "\t".join([str(x) for x in alt_pred]))

			ind1_pred = ind1_map()
			mse, spearman = comparePreds(np.array(alt_pred), np.array(ind1_pred))

			with open(f'IDW_modifying_variant_comparisons.txt', 'a') as comparison_out:
				comparison_out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr, variant_pos, window_start, individual, mse, spearman))
			
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

def alternate_map():
	args = parse_args()
	ref_fasta_open = pysam.Fastafile(f'{args.fasta}')
	ref_seq = ref_fasta_open.fetch(chr, int(window_start), int(window_start)+2**20).upper()

	pos_in_window = (int(variant_pos)-1)-int(window_start)
	alt_seq = ref_seq[0:pos_in_window] + alt_allele + ref_seq[pos_in_window+1:]
	alt_pred = runAkitaPreds(alt_seq)
	alt_pred_hff = alt_pred[:,:,0][0]
	return alt_pred_hff
	
def ind1_map():
	ind1_file = open(f'../predictions/HFF/3d_predictions_{individual}.txt', 'r')
	ind1_lines = ind1_file.readlines()
	for line in ind1_lines:
		if chr in line and window_start in line:
			ind1_line = line.split('\t')
			ind1_pred = list(map(float,ind1_line[2:]))
	return ind1_pred

def comparePreds(ind1, ind2):
		mse = np.mean(np.square(ind1 - ind2))
		spearman = stats.spearmanr(ind1, ind2)[0]
		return (mse, spearman)

if __name__ == '__main__':
	main()
