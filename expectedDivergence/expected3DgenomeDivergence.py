#!/usr/bin/env python

import os
import sys
import json
import configparser
import subprocess

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


def main():


    return

def setupAkitaModel():
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
    anc_indiv = 'hsmrca_ancestral'
    in_file_loc_anc_indiv = find_inFileLoc(anc_indiv, chrm)
    in_file_loc_mod_indiv = find_inFileLoc(mod_indiv, chrm)





    return

if __name__ == '__main__':
    main()