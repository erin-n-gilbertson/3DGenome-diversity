#!/usr/bin/env python

# Will calculate pairwise sequence similarities across the genome in sliding 1 Mb windows for pairs of individual
# Input: takes two individual ids  (example: python runSeqComparisons.py "vindija AFR_ESN_female_HG03105")
# Output: returns a file with sequence comparisons and coverage between the individuals

print("Started python",flush=True)
import numpy as np
import sys
import pysam
import configparser
import subprocess

configfile_name = sys.argv[2]

config = configparser.ConfigParser()
config.read(configfile_name)

### find file location for the individuals considered ###

def find_fastaFiles(indiv, chrm):

    fasta_dir = config["PATH"]["FASTA_PATH"]

    in_file_loc = '%s/%s/%s_%s_hg38_full.fa' % (fasta_dir,indiv,chrm,indiv.split("_")[3])

    indiv_fasta_open = pysam.Fastafile(in_file_loc)

    return indiv_fasta_open


#####################################

# input individuals
indivname1=sys.argv[1].strip().split(" ")[0]
indivname2=sys.argv[1].strip().split(" ")[1]

# read in chr/pos file over genome
GENOME_CHUNKS=config["FILE"]["GENOME_CHUNKS"]
chunks = {}
with open(GENOME_CHUNKS) as f:
    for line in f:
        chunk = line.strip().split("\t")
        chunks[chunk[0]] = [int(x) for x in chunk[1].split(",")]
f.close()

print("Indiv1 = %s, Indiv2 = %s" % (indivname1, indivname2),flush=True)

out_path = config["PATH"]["OUT_PATH"]
f_out = open("%s/SeqComps_%s_vs_%s.txt" % (out_path,indivname1,indivname2),"w")
f_out.write("%s\t%s\t%s\t%s\t%s\n" % ('chrm','start_loc',indivname1 + '_coverage',indivname2 + '_coverage', 'seqComp_raw'))

# Loop over chrms and positions
for chrm,pos_list in chunks.items():
    print("On chrom = %s" % chrm,flush=True)
    try:
        indiv1_fasta_open = find_fastaFiles(indivname1, chrm)
        indiv2_fasta_open = find_fastaFiles(indivname2, chrm)
    except OSError:
        print("Failed on chr: %s" % chrm)
        continue

    for start_loc in pos_list:
        print("starting sequence comparisons on %s" % start_loc,flush=True)
        try: # some input start locations won't work because when + 1Mb they are past the end of the chromosome stop
            # Fetch the fasta sequence
            indiv1_seq = indiv1_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()
            indiv2_seq = indiv2_fasta_open.fetch(chrm, start_loc, start_loc+2**20).upper()
            print("length of seq1: %s, length of seq2: %s" % (str(len(indiv1_seq)), str(len(indiv2_seq))))
            # calculate coverage
            indiv1_coverage = np.mean([ 0 if s == "N" else 1 for s in indiv1_seq])
            indiv2_coverage = np.mean([ 0 if s == "N" else 1 for s in indiv2_seq])
            # calculate sequence comparisons
            if len(indiv1_seq) !=0 :
                seqComp_raw = sum([1 if i1 == i2 else 0 for i1,i2 in zip(indiv1_seq,indiv2_seq)])/len(indiv1_seq)
            else:
                seqComp_raw = 'na'
            # write output to files
            f_out.write("%s\t%s\t%s\t%s\t%s\n" % (chrm,start_loc,indiv1_coverage,indiv2_coverage, seqComp_raw))
        except ValueError:
            print("FAILED: %s at %s" % (chrm, start_loc))
            continue

f_out.close()
