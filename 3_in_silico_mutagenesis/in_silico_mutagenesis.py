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


BASE_PATH = '/wynton/group/capra/projects/modern_human_3Dgenome'

BIN_PATH = os.path.join(BASE_PATH, "bin")  # where my scripts live
DATA_PATH = os.path.join(BASE_PATH, "data")  # where I dump new data 
RESULTS_PATH = os.path.join(BASE_PATH, "results")  # where I analyze results

FASTA_PATH = os.path.join(DATA_PATH, "genomes")
MODEL_DIR = '%s/basenji/manuscripts/akita' % BIN_PATH

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

		lines = input.readlines()
		if args.end is None:
			file_end = len(lines)
		elif args.end is not None:
			file_end = args.end
		print("args.start: %s" % args.start)	
		print("file_end: %s" % file_end)
		lines = lines[args.start:file_end]
		for line in lines:
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
			if alt_allele == "??":
				print("?? allele")
				df.append({'chr':chr,'pos':variant_pos,'window':int(window_start),'ref':ref_allele,'alt':alt_allele,'mse':'NA','1-pearson':'NA','1-spearman':'NA'})
			elif (len(alt_allele) > 1) or (len(ref_allele) > 1):
				print("indel allele")
				df.append({'chr':chr,'pos':variant_pos,'window':int(window_start),'ref':ref_allele,'alt':alt_allele,'mse':'NA','1-pearson':'NA','1-spearman':'NA'})
			else:
				print('regular allele')
				print("anc: %s, alt: %s" % (ref_allele, alt_allele))
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
		
	with open('%s/params.json' % MODEL_DIR) as params_file:
		params = json.load(params_file)
		params_model = params['model']
		params_train = params['train']

	global seqnn_model
	seqnn_model = seqnn.SeqNN(params_model)
	seqnn_model.restore('%s/model_best.h5' % MODEL_DIR)
	
def runAkitaPreds(seq):
	print(len(seq))
	print("2**20 = %s" % str(2**20))
	if len(seq) != 2**20: raise ValueError('len(seq) != seq_length')
	seq_1hot = dna_io.dna_1hot(seq)
	test_pred_from_seq = seqnn_model.model.predict(np.expand_dims(seq_1hot,0))
	return test_pred_from_seq

def find_inFileLoc(indiv, chrm):
    if indiv == 'human_archaic_ancestor':
        seq_dir = os.path.join(FASTA_PATH, "human_archaic_ancestor")
        file_path = os.path.join(seq_dir, ("human_archaic_ancestor_in_hg38_%s.fasta" % chrm))
    if indiv == 'hg38_reference':
        file_path = os.path.join(FASTA_PATH, ("hg38_reference.fa"))
    return file_path

def get_contact_maps():
	args = parse_args()
	in_file_loc_anc_indiv = find_inFileLoc(f'{args.fasta}', chr)
	ref_fasta_open = pysam.Fastafile(in_file_loc_anc_indiv)
	ref_seq = ref_fasta_open.fetch(chr, int(window_start), int(window_start)+2**20).upper()
	ref_pred = runAkitaPreds(ref_seq)
	ref_pred_hff = ref_pred[:,:,0][0]

	print("window start: %s" % window_start)
	print("variant position: %s" % variant_pos)
	print("length ref seq: %s" % len(ref_seq))

	pos_in_window = (int(variant_pos)-1)-int(window_start)
	if int(variant_pos) == int(window_start):
		print('variant pos == window start')
		alt_seq = alt_allele + ref_seq[1:]
	else:
		print('variant pos != window start')
		alt_seq = ref_seq[0:pos_in_window] + alt_allele + ref_seq[pos_in_window+1:]
	print("length alt seq: %s" % len(alt_seq))
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
