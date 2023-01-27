#!/usr/bin/env python

'''
Generates empirical distribution for the null hypothesis so we can calculate p values and enrichment statistics for observed counts
'''

import numpy as np
import pandas as pd
from pybedtools import BedTool
import sys
import argparse
from pandas.errors import EmptyDataError

ontology_results_dir = "../results/ontologyEnrichment"
ontology_data_dir = "../data/ontologyEnrichment"
tads = pd.read_csv(f"{ontology_data_dir}/4DN_HFF_MicroC_hg19_processedTADs.bed",sep="\t",header=None)
tads_bed = BedTool.from_dataframe(tads)
refseq = BedTool(f'{ontology_data_dir}/refSeqGeneCoordinatesSimplifiedProteinCoding.bed')

def main():
    args = parse_args()
    archaic_set = args.set
    vars_3d = getVariantsToShuffle(archaic_set)
    ontology = getGeneOntologyLinks(archaic_set)
    ontology_empiric = []
    while len(ontology_empiric) < args.number:
        empiric_counts = shuffleAndCalculateOntologyOverlap(vars_3d, ontology)
        if empiric_counts is not None:
            ontology_empiric.append(empiric_counts)
    if args.index is not None:
        outfilename = f"{args.outfile}{args.set}_{args.index}.tsv"
    else:
        outfilename = f"{args.outfile}{args.set}.tsv"
    pd.DataFrame(ontology_empiric).T.to_csv(outfilename, index=True, header=False, sep="\t") 

def parse_args():
    print("Reading in arguments...", flush=True)
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    parser.add_argument("-i", "--index", help="index for file name (if subsetting file for parallelization", type=str) 
    parser.add_argument("-n", "--number", help="number of empiric observations to generate", type=int, default=5000) 
    required.add_argument("-o", "--outfile", help="path/name for outfile to write to", type=str, required=True) 
    required.add_argument("-s", "--set", help="set name for 3d modifying variants to shuffle (options = intersectNeanderthal, unionNeanderthal,denisova", type=str, required=True) 
    return parser.parse_args()

def getVariantsToShuffle(archaic_set):
    allTopVars = pd.read_csv("../results/AH-MH_divergedWindows_inSilicoMut.tsv",sep="\t")
    allTopVars['AH'] = [set([y.strip().split(":")[0] for y in x.strip().split(";")]) for x in allTopVars['AH_inSilicoMut'].values]
    allTopVars['intersectNeanderthal'] = [('V' in x) and ('A' in x) and ('C' in x) for x in allTopVars['AH']]
    allTopVars['unionNeanderthal'] = [('V' in x) or ('A' in x) or ('C' in x) for x in allTopVars['AH']]
    allTopVars['denisova'] = [('D' in x) for x in allTopVars['AH']]
    vars_3d = allTopVars[allTopVars[archaic_set]]
    vars_3d['start'] = vars_3d['pos'] - 1
    vars_3d = BedTool().from_dataframe(vars_3d[['chr','start','pos']]).sort()
    return vars_3d

def getGeneOntologyLinks(archaic_set):
    hpo_observed_ontology_df = pd.read_csv(f'{ontology_results_dir}/observedCounts/hpo_ontology_counts_{archaic_set}.tsv.gz',sep="\t",index_col=0)
    gwas_observed_ontology_df = pd.read_csv(f'{ontology_results_dir}/observedCounts/gwas_ontology_counts_{archaic_set}.tsv.gz',sep="\t",index_col=0)
    observed_ontology_df = pd.concat([hpo_observed_ontology_df, gwas_observed_ontology_df])
    terms = list(observed_ontology_df.index[observed_ontology_df.sum(axis=1) > 0])
    ontology = {}
    for o in ['GWAS_Catalog_2019','Human_Phenotype_Ontology']:
        file = open(f'{ontology_data_dir}/{o}.txt', 'r')
        Lines = file.readlines()
        for line in Lines:
            line = line.strip().split("\t")
            if o + ": " + line[0] in terms:
                ontology[o + ": " + line[0]] = line[2:]
        file.close()
    return ontology

def shuffleAndCalculateOntologyOverlap(vars_3d, ontology):
    try:
        vars_3d_shuffle = vars_3d.shuffle(genome='hg19',incl = "../data/regionsWithFullCoverage.bed")
        intersect = tads_bed.intersect(vars_3d_shuffle,wo=True).intersect(refseq,loj=True).to_dataframe(disable_auto_names=True, header=None)
        genes = list([x for x in intersect[10] if str(x) != '.'])
        empiric_counts = {}
        for i,r in ontology.items():
            empiric_counts[i]  = 0
            for g in genes:
                if g in r:
                    empiric_counts[i]+=1
        return empiric_counts
    except EmptyDataError:
        return None

if __name__ == '__main__':
    main()
