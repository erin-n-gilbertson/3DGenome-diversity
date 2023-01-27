#!/usr/bin/env python
'''
Conducts in silico mutagenesis over a file of 1 Mb windows across the genome

usage: inSilicoMutagenesis.py [-h] -f FILE [-s STARTLINE] [-e ENDLINE] -o
                              OUTFOLDER
optional arguments:  -h, --help            show this help message and exit
  -s STARTLINE, --startline STARTLINE
                        start line in bedfile if you want to use a subset of
                        lines (Default = 0)
  -e ENDLINE, --endline ENDLINE
                          end line in bedfile if you want to use a subset of
                          lines (Default = last line in file)required arguments:
  -f FILE, --file FILE  path to file
  -o OUTFOLDER, --outfolder OUTFOLDER
                          path of folder to write to
'''
print("Reading in dependencies...",flush=True)
import argparse
import random 
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import ast
import pysam

def main():
    '''
    Parses arguments, reads in file, loads Akita, loops over infile entry to:
    make seqs with and without in silico mutagenesis, run Akita predictions, 
    compare predictions, write file.
    '''
    args = parse_args()
    df_pos = readInFile(args.file, args.startline, args.endline)
    loadAkita()
    for _,r in df_pos.iterrows(): # For each file entry
        for ah in r['AH']: #[x.split("_")[0] for x in r['AH']]:
            print(f"On {r['chr']}:{r['windowStartPos']} ({ah})...",flush=True)
            insilicoMutagenesis(r['chr'], r['windowStartPos'], ah, args.outfolder)


def parse_args():
    print("Reading in arguments...", flush=True)
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--file", help="path to file", required=True)
    parser.add_argument("-s", "--startline", help="start line in bedfile if you want to use a subset of lines (Default = 0)", default = 0, type=int) 
    parser.add_argument("-e", "--endline", help="end line in bedfile if you want to use a subset of lines (Default = last line in file)", type=int) 
    required.add_argument("-o", "--outfolder", help="path of folder to write to", type=str, required=True) 
    return parser.parse_args()


def insilicoMutagenesis(chrm, start_loc, arc_indiv, outfolder):
    afr_20 = ['AFR_ESN_female_HG03105',
        'AFR_GWD_female_HG03539',
        'AFR_LWK_female_NA19378',
        'AFR_MSL_female_HG03212',
        'AFR_ESN_female_HG03499',
        'AFR_ESN_female_HG03511',
        'AFR_ESN_female_HG03514',
        'AFR_GWD_female_HG03025',
        'AFR_GWD_female_HG03028',
        'AFR_LWK_female_NA19017',
        'AFR_LWK_female_NA19434',
        'AFR_LWK_female_NA19445',
        'AFR_GWD_female_HG03040',
        'AFR_ESN_female_HG02922',
        'AFR_MSL_female_HG03086',
        'AFR_MSL_female_HG03085',
        'AFR_GWD_female_HG03046',
        'AFR_MSL_female_HG03437',
        'AFR_MSL_female_HG03378',
        'AFR_LWK_female_NA19019',]

    in_file_loc_arc_indiv = find_inFileLoc(arc_indiv, chrm)
    arc_indiv_fasta_open = pysam.Fastafile(in_file_loc_arc_indiv)
    arc_indiv_seq = arc_indiv_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()

    mask_fasta_open = pysam.Fastafile('/gpfs51/dors2/capra_lab/users/rinkerd/projects/3DNeand/data/genomes/masked_hg19_reference/%s_hg19_archaic.masked.fa' % chrm) #for the masked
    masked_seq = mask_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()

    human19_fasta_open = pysam.Fastafile('/dors/capra_lab/data/dna/human/hg19/%s.fa' % chrm)
    human19_seq = human19_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()

    arc_masked = "".join(["N" if m == "N" else "N" if ref == "N" else arc for m, ref, arc in zip(masked_seq, human19_seq, arc_indiv_seq)])
    afr_seq = {}
    for afr_indiv in afr_20:
        in_file_loc_afr_indiv = find_inFileLoc(afr_indiv, chrm)
        afr_indiv_fasta_open = pysam.Fastafile(in_file_loc_afr_indiv)
        afr_indiv_seq = afr_indiv_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()
        afr_masked = "".join(["N" if m == "N" else "N" if ref == "N" else "N" if arc == "N" else afr for m, ref, arc, afr in zip(masked_seq, human19_seq, arc_indiv_seq, afr_indiv_seq)])
        afr_seq[afr_indiv] = afr_masked
    diff = [all([v[i]!=arc_masked[i] for _,v in afr_seq.items()]) for i in range(2**20)]
    diff_num = sum(diff)
    print("Number of ISM positions: " + str(diff_num))
    afr_filled = "".join([r if m == "N" else r if s == "N" else s for r, m, s in zip(human19_seq, masked_seq, afr_seq['AFR_ESN_female_HG03105'])])
    afr_pred  = runAkitaPreds(afr_filled)

    df = []
    for i in range(len(afr_masked)):
        if diff[i]:
            print(str(i))
            afr_withArcMut_filled = afr_filled[:i] + arc_indiv_seq[i] + afr_filled[i+1:]
            afr_withArcMut_pred = runAkitaPreds(afr_withArcMut_filled)
            spearman, pearson, mse = compare3Dmaps(afr_pred, afr_withArcMut_pred)
            df.append({'chr':chrm, 'pos':start_loc+i+1,'1-spearman':spearman,'1-pearson':pearson,'mse':mse,'afr_allele':afr_indiv_seq[i], 'arc_allele':arc_indiv_seq[i] })
    df = pd.DataFrame(df)[['chr','pos','afr_allele','arc_allele','1-spearman','1-pearson','mse']]
    df.to_csv(f"{outfolder}/inSilicoMut_%s_%s_%s.tsv" % (arc_indiv, chrm, start_loc),sep="\t",index=False)

def runAkitaPreds(seq):
    ''' Runs Akita prediction on `seq` (ACTG string of 2^20 length) and outputs flattened vector for 2D map predictions (length = 99,681) '''
    if len(seq) != 2**20: raise ValueError('len(seq) != seq_length')
    seq_1hot = one_hot_encode(seq) #seq_1hot = dna_io.dna_1hot(seq)
    test_pred_from_seq = seqnn_model.model.predict(np.expand_dims(seq_1hot,0))
    return test_pred_from_seq[:,:,0][0] # Returns the HFF cell type, use [:,:,1][0] for ESC

def readInFile(infile, startline, endline):
    ''' Reads the divergent Windows file into a data frame and keeps only the lines between [startline, endline)
        Parameters:
            infile (str): path to file of positions
            Looks like this:
                chr     windowStartPos  AH
                chr1    38273024        ['denisova']
                chr1    56623104        ['denisova']
                chr1    66584576        ['vindija', 'chagyrskaya', 'altai', 'denisova']
            startline (int): start line in file
            endline (int or None): end line in file (if None will use last line)
        Returns:
            df_pos (df): data frame with lines from file specified
    '''
    print('Reading in  input file...',flush=True)
    df = pd.read_csv(infile,sep="\t")
    df['AH'] = [[x.strip() for x in ast.literal_eval(r)] for r in df['AH'].values]
    if endline is None:
        return df.iloc[startline:]
    else:
        return df.iloc[startline:endline]

def loadAkita():
    '''Read in dependencies, hg38 and Akita model (as global vars)'''
    print("Loading Akita...",flush=True)
    os.environ["KMP_WARNINGS"] = "FALSE" 
    import json
    os.environ["CUDA_VISIBLE_DEVICES"] = '-1' ### run on CPU
    import tensorflow as tf
    if tf.__version__[0] == '1':
        tf.compat.v1.enable_eager_execution()
    from basenji import dataset, seqnn#, dna_io
    ### load params, specify model ###
    model_dir = './akita_model/'
    params_file = model_dir+'params.json'
    model_file  = model_dir+'model_best.h5'
    with open(params_file) as params_open:
        params = json.load(params_open)
        params_model = params['model']
        params_train = params['train']
    global seqnn_model
    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore(model_file)

def find_inFileLoc(indiv, chrm):
    if (indiv == "vindija") | (indiv == "altai") | (indiv == "denisova") | (indiv == "chagyrskaya"):   # if archaic
        in_file_loc = '/gpfs51/dors2/capra_lab/users/rinkerd/projects/3DNeand/data/genomes/%s/%s%s_hg19_masked.fa' % (indiv,chrm,indiv)
    else:        
        in_file_loc = '/gpfs51/dors2/capra_lab/users/rinkerd/projects/3DNeand/data/genomes/1KG.individuals/%s/%s_%s_hg19_full.fa' % (indiv,chrm,indiv.split("_")[3])
    return in_file_loc

def one_hot_encode(sequence: str,
                alphabet: str = 'ACGT',
                neutral_alphabet: str = 'N',
                neutral_value = 0,
                dtype=np.float32) -> np.ndarray:
    '''Faster One-hot encode sequence. From Enformer paper and L. Gunsalus'''
    def to_uint8(string):
        return np.frombuffer(string.encode('ascii'), dtype=np.uint8)
    hash_table = np.zeros((np.iinfo(np.uint8).max, len(alphabet)), dtype=dtype)
    hash_table[to_uint8(alphabet)] = np.eye(len(alphabet), dtype=dtype)
    hash_table[to_uint8(neutral_alphabet)] = neutral_value
    hash_table = hash_table.astype(dtype)
    return hash_table[to_uint8(sequence)]

def compare3Dmaps(pred1, pred2):
    ''' 
    Compares the contact maps (from the output flattened vectors in runAkitaPreds) to each other
    Parameters:
        pred1 (vector): contact map flattened vector #1
        pred2 (vector): contact map flattened vector #2
    Returns: 
        Tuple with n=3 comparisons (all scaled to values where higher = more difference between maps)
        (1 - spearman,  1 - pearson,  MSE)
    '''
    spearman = 1 - stats.spearmanr(pred1 ,pred2 )[0]
    pearson = 1 - stats.pearsonr(pred1 ,pred2 )[0]
    mse = np.mean(np.square(pred1  - pred2 ))
    return spearman, pearson, mse

if __name__ == '__main__':
    main()
