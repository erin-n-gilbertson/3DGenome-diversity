#!/usr/bin/env python

'''
Used to generate an expected number of false positives for ontology enrichment in order to calculate FDR-corrected P values (q)
'''


import pandas as pd
import numpy as np
import argparse
import gzip

def main():
    args = parse_args()
    f = gzip.open(args.file, 'r')
    lines = f.readlines()
    out_df = []
    for count, line in enumerate(lines):
        if count != args.index:
            continue
        line = np.array(line.decode().strip().split("\t"))
        label = line[0]
        if args.number is not None:
            line = line[1:args.number+1]
        else:
            line = line[1:]
        pvals = []
        for i, v in enumerate(line):
            if i % 500 == 0:
                print(i,flush=True)
            p = sum(line[:i] >= v) + sum(line[i+1:] >= v)
            p = (p+1)/len(line)
            pvals.append(p)
        out_df.append({**{'label': label},**dict(zip(range(len(pvals)),pvals))})
    out_df = pd.DataFrame(out_df)
    out_df.to_csv(f"{args.outfile}_{args.index}.tsv",sep="\t",header=None, index=None)

def parse_args():
    print("Reading in arguments...", flush=True)
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--index", help="index for which row in input file to consider", type=int,required=True) 
    required.add_argument("-o", "--outfile", help="path/name for outfile to write to", type=str, required=True) 
    required.add_argument("-f", "--file", help="input file for empirical counts", type=str, required=True) 
    parser.add_argument("-n", "--number", help="number of empirical counts/observations per row to consider. This gets expensive as you increase!!", type=int) 
    return parser.parse_args()


if __name__ == '__main__':
    main()

