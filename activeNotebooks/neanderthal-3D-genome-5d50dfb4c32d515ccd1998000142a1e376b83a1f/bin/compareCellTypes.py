#!/usr/bin/env python
""" Runs Akita on positions across the genome to compare HFF and ESC results. 
    Example usage: python compareCellTypes.py --i1 vindija --i2 AFR_ESN_female_HG03105 --dir cellTypeComps_intermediate --pos_list ../data/regionsWithFullCoverage.csv --rows 1-5
                   python compareCellTypes.py --i1 vindija --i2 AFR_ESN_female_HG03105 --dir cellTypeComps_intermediate --pos_list ../data/regionsWithFullCoverage.csv --rows 6-10
                   ... (can parallelize here)
                   ls cellTypeComps_intermediate/* > filesToJoin.txt
                   python compareCellTypes.py --combine filesToJoin.txt --out 3dcompAcrossCellTypes_AFR_ESN_female_HG03105_vs_vindija
                   gzip 3dcompAcrossCellTypes_AFR_ESN_female_HG03105_vs_vindija.tsv
                   rm -rf cellTypeComps_intermediate/ filesToJoin.txt
"""
print("Loading dependencies... Tensorflow version:",flush=True)
import os
import sys
import argparse

###### Import to run Akita ######
os.environ["KMP_WARNINGS"] = "FALSE" 
import numpy as np
np.random.seed(42)
import pandas as pd
import json
os.environ["CUDA_VISIBLE_DEVICES"] = '-1' ### run on CPU
import tensorflow as tf
print(tf.__version__,flush=True)
if tf.__version__[0] == '1':
    tf.compat.v1.enable_eager_execution()
import pysam
import scipy.stats as stats
from basenji import dataset, dna_io, seqnn


def main(arguments):
    print("Reading in args...",flush=True)
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--i1", help="Name of Individual #1", type=str, default=argparse.SUPPRESS)
    parser.add_argument("--i2", help="Name of Individual #2", type=str, default=argparse.SUPPRESS)
    parser.add_argument("--dir", help="Directory name", type=str, default="default_dir")
    parser.add_argument('--pos_list', help="File with positions to run Akita on (comma sep file with chr,pos_start in each row)", type=str,default=False)
    parser.add_argument('--rows', help="Line numbers (0 indexed) of rows in pos_list file to run (helpful for parallelization), format '#-#', Default: uses all lines in file", type=str, default=False)
    parser.add_argument('--combine', help="List of files to combine", type=str,default=False)
    parser.add_argument('--out', help="Name of output file", type=str, default="compCellTypes")
    args = parser.parse_args(arguments)
    print(args,flush=True)
    if args.combine == False:
        create_dir(args.dir)
        df_pos = readInFileWithPositions(args.pos_list, args.rows)
        runAkita(args.i1, args.i2, args.dir,df_pos,args.rows, args.out)
    else:
        combineFiles(args.combine, args.out)
######################################################################

def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
        print("Created Directory : ", dir)

def readInFileWithPositions(pos_list, rows):
    print('Reading in file of positions...',flush=True)
    df = pd.read_csv(pos_list,header=None)
    row_start = int(rows.split("-")[0])
    row_end = int(rows.split("-")[1])
    df = df.iloc[row_start:row_end+1]
    return df.sort_values([0,1])

def runAkita(i1, i2, dir, df_pos, rows, out_str):
    print("Loading Akita...",flush=True)
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

    ### restore model ###
    seqnn_model.restore(model_file)

    ### read data parameters ###
    data_stats_file = model_dir+'/statistics.json'
    with open(data_stats_file) as data_stats_open:
        data_stats = json.load(data_stats_open)
    seq_length = data_stats['seq_length']
    target_length = data_stats['target_length']
    hic_diags =  data_stats['diagonal_offset']
    target_crop = data_stats['crop_bp'] // data_stats['pool_width']
    target_length1 = data_stats['seq_length'] // data_stats['pool_width']

    ### go through positions in df_pos list
    divergence_df = []
    for _,r in df_pos.iterrows():
        print(r[0],r[1],flush=True)
        indiv1 = extractFastaRunAkita(i1,r[0],int(r[1]))
        indiv2 = extractFastaRunAkita(i2,r[0],int(r[1]))
        hff_spearman = 1 - stats.spearmanr(indiv1.pred[:,:,0][0] ,indiv2.pred[:,:,0][0] )[0]
        hff_pearson = 1 - stats.pearsonr(indiv1.pred[:,:,0][0] ,indiv2.pred[:,:,0][0] )[0]
        hff_mse = np.mean(np.square(indiv1.pred[:,:,0][0]  - indiv2.pred[:,:,0][0] ))
        esc_spearman = 1 - stats.spearmanr(indiv1.pred[:,:,1][0] ,indiv2.pred[:,:,1][0])[0]
        esc_pearson = 1 - stats.pearsonr(indiv1.pred[:,:,1][0],indiv2.pred[:,:,1][0])[0]
        esc_mse = np.mean(np.square(indiv1.pred[:,:,1][0] - indiv2.pred[:,:,1][0]))
        divergence_df.append([r[0],r[1],hff_spearman, hff_pearson, hff_mse,esc_spearman,esc_pearson, esc_mse])
    divergence_df = pd.DataFrame(divergence_df)
    divergence_df.columns = ['chr','windowStartPos','hff_spearman', 'hff_pearson', 'hff_mse','esc_spearman','esc_pearson', 'esc_mse']
    divergence_df = divergence_df.sort_values(['chr','windowStartPos'])
    divergence_df.to_csv(f"{dir}/{out_str}_{rows}.tsv",sep="\t",index=False)
    print("Done",flush=True)
    
def find_inFileLoc(indiv, chrm):
    if (indiv == "vindija") | (indiv == "altai") | (indiv == "denisova") | (indiv == "chagyrskaya"):   # if archaic
        in_file_loc = '/gpfs51/dors2/capra_lab/users/rinkerd/projects/3DNeand/data/genomes/%s/%s%s_hg19_masked.fa' % (indiv,chrm,indiv)
    elif 'SAS' not in indiv:
        in_file_loc = '/gpfs51/dors2/capra_lab/users/rinkerd/projects/3DNeand/data/genomes/1KG.individuals/%s/%s_%s_hg19_full.fa' % (indiv,chrm,indiv.split("_")[3])
    else:
        in_file_loc = '/dors/capra_lab/users/erin/RotationProject_Akita/data/genomes/1KG/%s/%s_%s_hg19_full.fa' % (indiv,chrm,indiv.split("_")[3])

    return in_file_loc

### Run Akita 3d genome predictions ###
def runAkitaPreds(seq):
    if len(seq) != 2**20: raise ValueError('len(seq) != seq_length')
    seq_1hot = dna_io.dna_1hot(seq)
    test_pred_from_seq = seqnn_model.model.predict(np.expand_dims(seq_1hot,0))
    return test_pred_from_seq

### Object that gets the sequence and 3d genome predictions for a window for a given individual ###
class extractFastaRunAkita:
    def __init__(self, indiv, chrm, start_loc):
        in_file_loc_indiv = find_inFileLoc(indiv, chrm)
        try:
            indiv_fasta_open = pysam.Fastafile(in_file_loc_indiv)
        except OSError:
            in_file_loc_indiv = '/dors/capra_lab/users/erin/RotationProject_Akita/data/genomes/1KG/%s/%s_%s_hg19_full.fa' % (indiv,chrm,indiv.split("_")[3])
            indiv_fasta_open = pysam.Fastafile(in_file_loc_indiv)

        mask_fasta_open = pysam.Fastafile('/gpfs51/dors2/capra_lab/users/rinkerd/projects/3DNeand/data/genomes/masked_hg19_reference/%s_hg19_archaic.masked.fa' % chrm) #for the masked
        human19_fasta_open = pysam.Fastafile('/dors/capra_lab/data/dna/human/hg19/%s.fa' % chrm)

        #extract sequences
        indiv_seq = indiv_fasta_open.fetch(chrm, int(start_loc), start_loc+2**20).upper()
        masked_seq = mask_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper() #for the masked 
        human19_seq = human19_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()

        # important harmonization step! then run predictions on the harmonized sequence
        indiv_fillMissing_seq = "".join([r if m == "N" else r if s == "N" else s for r, m, s in zip(human19_seq, masked_seq, indiv_seq)])
        self.seq = indiv_fillMissing_seq
        indiv_pred  = runAkitaPreds(indiv_fillMissing_seq)
        self.pred = indiv_pred # save all cell types!
        
def combineFiles(combine_file, out_str):
    print("Combining Files...",flush=True)
    f = open(combine_file, 'r')
    lines = f.readlines()
    df = pd.DataFrame(columns=['chr','windowStartPos','hff_spearman', 'hff_pearson', 'hff_mse','esc_spearman','esc_pearson', 'esc_mse'])
    for l in lines:
        df = pd.concat([df,pd.read_csv(l.strip(),sep="\t")])
    df.sort_values(['chr','windowStartPos']).to_csv(out_str + ".tsv",sep="\t",index=False)
    f.close()
    print('Done!',flush=True)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))