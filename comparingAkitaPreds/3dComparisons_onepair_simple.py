#!/usr/bin/env python

print("Started python",flush=True)
import numpy as np
from scipy import stats
#from cooltools.lib.numutils import set_diag
import sys
from itertools import chain
import warnings
import os.path
from os import path
import gzip

configfile_name = sys.argv[2]

config = configparser.ConfigParser()
config.read(configfile_name)

### for converting from flattened upper-triangluar vector to symmetric matrix  ###

# def from_upper_triu(vector_repr, matrix_len, num_diags):
    # z = np.zeros((matrix_len,matrix_len))
    # triu_tup = np.triu_indices(matrix_len,num_diags)
    # z[triu_tup] = vector_repr
    # for i in range(-num_diags+1,num_diags):
        # set_diag(z, np.nan, i)
    # return z + z.T



def comparePreds(indiv1, indiv2):
    mse = np.mean(np.square(indiv1 - indiv2))
    spearman = stats.spearmanr(indiv1, indiv2)[0]
    return (mse, spearman)

#####################################

# Useful for slurm scripts
indivname1=sys.argv[1].split(" ")[0]
indivname2=sys.argv[1].split(" ")[1]

#indivname1=sys.argv[1]
#indivname2=sys.argv[2]
#path2indiv1=sys.argv[3]
#path2indiv2=sys.argv[4]



print("Indiv1 = %s, Indiv2 = %s" % (indivname1, indivname2),flush=True)


in_file_loc = config["PATH"]["INPUT_PATH"]


# if using gzip file read binary and use .decode() for l1 & l2 below (4x)

f1 = open("%s/3dpreds_%s.txt" % (in_file_loc, indivname1),"rb")
f2 = open("%s/3dpreds_%s.txt" % (in_file_loc, indivname2),"rb")

f_out = open("%s/3dcomp_%s_vs_%s.txt" % (config["PATH"]["OUT_PATH"],indivname1,indivname2),"w")

# f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr", "windowStartPos", "mse", "spearman", "triangle_mse","triangle_spearman", "insulation_spearman"))
f_out.write("%s\t%s\t%s\t%s\n" % ("chr", "windowStartPos", "mse", "spearman"))

lines1 = f1.readlines()
lines2 = f2.readlines()
shift1 = 0
shift2 = 0
print(len(lines1))
print(len(lines2))
for i in range(max(len(lines1),len(lines2))): # loop through indexes of the longest file

  print(i, flush=True)

  # Handle instances when one file has extra regions at the end
  if (i + shift1) >= len(lines1):
    l2 = lines2[i + shift2].decode()
    indiv2 = l2.strip().split("\t")
    indiv2_chr = indiv2[0]
    indiv2_pos = indiv2[1]
    print("File 1 missing %s:%s" % (indiv2_chr, indiv2_pos))
    continue
  if (i + shift2) >= len(lines2):
    l1 = lines1[i + shift1].decode()
    indiv1 = l1.strip().split("\t")
    indiv1_chr = indiv1[0]
    indiv1_pos = indiv1[1]
    print("File 2 missing %s:%s" % (indiv1_chr, indiv1_pos))
    continue

  # Read in files and separate out their chrm/position/3dpred info
  l1 = lines1[i + shift1].decode()
  l2 = lines2[i + shift2].decode()
  indiv1 = l1.strip().split("\t")
  indiv1_chr = indiv1[0]
  indiv1_pos = indiv1[1]
  indiv1 = list(map(float,indiv1[2:]))
  indiv2 = l2.strip().split("\t")
  indiv2_chr = indiv2[0]
  indiv2_pos = indiv2[1]
  indiv2 = list(map(float,indiv2[2:]))

  # If the position at the current index is the same in both files, just simply output the comparison data
  if (indiv1_chr == indiv2_chr) and (indiv1_pos == indiv2_pos):
    mse, spearman = comparePreds(np.array(indiv1), np.array(indiv2))
    f_out.write("%s\t%s\t%s\t%s\n" % (indiv1_chr, indiv1_pos, mse, spearman))
        # mse, spearman, triangle_mse,triangle_spearman, insulation_spearman = comparePreds(np.array(indiv1), np.array(indiv2))
    # f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (indiv1_chr, indiv1_pos, mse, spearman, triangle_mse, triangle_spearman, insulation_spearman))

  # If the position is not the same in both files, traverse one file by decrementing the shift variable
  elif (indiv1_chr == indiv2_chr):
    if int(indiv1_pos) > int(indiv2_pos):
      shift1-=1
      print("File 1 missing %s:%s" % (indiv2_chr, indiv2_pos))
    else:
      shift2-=1
      print("File 2 missing %s:%s" % (indiv1_chr, indiv1_pos))
  else:
    indiv1_chr_sub = indiv1_chr.split("chr")[1]
    indiv2_chr_sub = indiv2_chr.split("chr")[1]
    # handling the X chromosome comparisons
    if indiv1_chr_sub == "X":
        indiv1_chr_sub = 23
    if indiv2_chr_sub == "X":
        indiv2_chr_sub = 23
    if int(indiv1_chr_sub) > int(indiv2_chr_sub):
      shift1-=1
      print("File 1 missing %s:%s" % (indiv2_chr, indiv2_pos))
    else:
      shift2-=1
      print("File 2 missing %s:%s" % (indiv1_chr, indiv1_pos))

f_out.close()
f1.close()
f2.close()
