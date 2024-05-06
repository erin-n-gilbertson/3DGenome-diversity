import pandas as pd
import pybedtools
import pysam

import argparse
import os
import json
from basenji import dataset, dna_io, seqnn
import scipy.stats as stats


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
        "--indiv", type=str, required=True,
        help="1KG individual to use")

    parser.add_argument(
        "--window_size_exponent", type=str, required=True,
        help="window size for comparison, represented as the exponent for 2**n that will be used to define window size and related characteristics")
    args = parser.parse_args()
    return args


def get_anc_seq(chr, window_start):
    anc_fasta_open = pysam.Fastafile('%s/genomes/human_archaic_ancestor/human_archaic_ancestor_in_hg38_%s.fasta' % (DATA_PATH, chr) )

    anc_seq = anc_fasta_open.fetch(chr, int(window_start), int(window_start)+2**20).upper()

    return anc_seq

def get_ref_seq_hg38(chr, window_start):
    ref_fasta_open = pysam.Fastafile('%s/genomes/hg38_reference.fa' % DATA_PATH)

    ref_seq = ref_fasta_open.fetch(chr, int(window_start), int(window_start)+2**20).upper()

    return ref_seq

def comparePreds(pred1, pred2):
    mse = np.mean(np.square(pred1 - pred2))
    spearman = stats.spearmanr(pred1, pred2)[0]
    divergence = 1 - spearman
    return (mse, divergence)

def from_upper_triu(vector_repr, matrix_len, num_diags):
    z = np.zeros((matrix_len,matrix_len))
    triu_tup = np.triu_indices(matrix_len,num_diags)
    z[triu_tup] = vector_repr
    for i in range(-num_diags+1,num_diags):
        set_diag(z, np.nan, i)
    return z + z.T

def sym_mat(pred, s=448):
    ind1_mat = from_upper_triu(pred, s, 2)
    # mask =  np.tri(ind1_mat.shape[0], k = -1) # bottom-half
    # ind1_mat = np.ma.array(ind1_mat, mask = mask).T # transpose

    return ind1_mat

def flatten(arry, diag_offset=2):
    flat = np.array([])
    rnum = diag_offset
    for row in arry:
        flat = np.concatenate([flat, row[rnum:]])
        rnum+=1
    return flat

def loadAkita():
    os.environ["CUDA_VISIBLE_DEVICES"] = '-1'

    import tensorflow as tf
    if tf.__version__[0] == '1':
        tf.compat.v1.enable_eager_execution()

    with open('%s/basenji/manuscripts/akita/params.json' % BIN_PATH) as params_file:
        params = json.load(params_file)
        params_model = params['model']
        params_train = params['train']
        
    global seqnn_model
    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore('%s/basenji/manuscripts/akita/model_best.h5' % BIN_PATH)

def runAkitaPreds(seq):
    if len(seq) != 2**20: raise ValueError('len(seq) != seq_length')
    seq_1hot = dna_io.dna_1hot(seq)
    test_pred_from_seq = seqnn_model.model.predict(np.expand_dims(seq_1hot,0))
    return test_pred_from_seq

# this function will load up a predicted map that I've already generated for the specified window. 
def load_1kg_map(chr, start, ind):
    with open('/wynton/group/capra/projects/modern_human_3Dgenome/data/akitaPreds/3dpreds/HFF_original/3dpreds_%s.txt' % (ind)) as file:
        lines = [ line.strip() for line in file ]
        for line in lines:
            if '%s' % (chr) in line and '%s' % (start) in line:
                ind2_vector = line.split('\t')
                ind2_vector = ind2_vector[2:]
    return ind2_vector



def main():
    args = parse_args()
    loadAkita()
    i = args.indiv
    e = int(args.window_size_exponent)
    window_size = 2**e
    mat_size = window_size/2048
    crop_size = mat_size - 2**(e-14)
    divisor = 2**(20-e) * 2
    num_subs = divisor - 1

    windows = pd.read_table('%s/intermediates/windows_to_keep.csv' % DATA_PATH, sep=',', index_col=0)
    
    df = []
    df.append(['window_start', 'window_end', 'mse', 'div'])
    for sw in range(num_subs):
        sub_start = 'sub_start_%s' % sw
        sub_end = 'sub_end_%s' % sw
        sub_mse = 'sub_mse_%s' % sw
        sub_div = 'sub_div_%s' % sw
        df.append([sub_start, sub_end, sub_mse, sub_div])


    for w in windows.index:
        wlist = []
        chr = windows.loc[w].chr
        window_start = windows.loc[w].windowStartPos
        print(chr,window_start)

        i_fasta = pysam.Fastafile('%s/genomes/1KG/%s/%s/%s_%s_hg38_full.fa' % (DATA_PATH, i.split('_')[0], i, chr, i.split('_')[-1]))
        
        anc_seq = get_anc_seq(chr, window_start)
        anc_pred = runAkitaPreds(anc_seq)
        anc_pred = anc_pred[:,:,0][0]

        sym_pred_anc = sym_mat(anc_pred)

        i_seq = i_fasta.fetch(chr, window_start, window_start+2**20).upper()
        i_pred = runAkitaPreds(i_seq)
        i_pred = i_pred[:,:,0][0]

        sym_pred_i = sym_mat(i_pred)
        mse, div = comparePreds(flatten(sym_pred_anc).astype('float32'), flatten(sym_pred_i).astype('float32'))

        wlist += [window_start, window_start + 2**20, mse, div]
        idx0, idx1 = 0, crop_size
        seq_idx0, seq_idx1 = (mat_size-crop_size)/2, (mat_size-crop_size)/2+crop_size
        for sw in range(num_subs):

            i_pred_sw = i_pred_sw[idx0:idx1, idx0:idx1]
            anc_pred_sw = anc_pred_sw[idx0:idx1, idx0:idx1]

            mse_sw, div_sw = comparePreds(flatten(anc_pred_sw).astype('float32'), flatten(i_pred_sw).astype('float32'))
            wlist += [seq_idx0*2048, seq_idx1*2048, mse_sw, div_sw]

            idx0 += crop_size/2
            idx1 += crop_size/2
            seq_idx0 += crop_size/2
            seq_idx1 += crop_size/2

        df.append(wlist)

    i_fasta.close()
    df = pd.DataFrame(df)
    df.to_csv(args.out, sep='\t', index=False)
    



if __name__ == '__main__':
    main()