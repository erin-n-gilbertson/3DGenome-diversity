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

### for converting from flattened upper-triangluar vector to symmetric matrix  ###

# def from_upper_triu(vector_repr, matrix_len, num_diags):
    # z = np.zeros((matrix_len,matrix_len))
    # triu_tup = np.triu_indices(matrix_len,num_diags)
    # z[triu_tup] = vector_repr
    # for i in range(-num_diags+1,num_diags):
        # set_diag(z, np.nan, i)
    # return z + z.T

### Insulation indicies ###

# def triangleScore(input_tuple):
    # mat = from_upper_triu(range(99681), 448, 2)
    # shift_i, triangle_i = input_tuple
    # return list(mat[tuple(i+shift_i for i in np.triu_indices(triangle_i,2))].astype(int))

# def triangleIndicies():
    # insulation_tri_indices = []

    # shift_max = 350 #448 - 3# target_length1_cropped - 3 #7-3
    # triangle_i = 4
    # while shift_max >= 0:
        # shift_i = 0
        # while shift_i < shift_max:
            # insulation_tri_indices.append((shift_i, triangle_i))
            # shift_i+=4 #!change to 1
        # shift_max-= 4 #!change to 1
        # triangle_i+= 4 #!change to 1
    # listOfIndexes = []
    # indexes= np.array(range(99681))
    # for i in insulation_tri_indices:
        # listOfIndexes.append(list(indexes[triangleScore(i)]))
    # return listOfIndexes

# print("creating insulation indicies",flush=True)
# listOfIndexes = triangleIndicies()
# print("done",flush=True)

# def get_insulation_track(mat, center_size, max_dist):
    # insulation_track = []
    # for b in range(0,len(mat),4): # change to no skipped spaces
        # if b <= center_size-1 or b >= len(mat)-center_size:
            # insulation_track.append(np.nan)
            # continue
        # with warnings.catch_warnings():
            # warnings.filterwarnings(action='ignore', message='Mean of empty slice')
            # focal_start = b-center_size-1
            # focal_end = b+center_size
            # center = mat[focal_start:focal_end, focal_start:focal_end]
            # center_mean = np.nanmean([np.exp(i) for i in list(chain(*center))])
            # upstream_b = focal_start-max_dist
            # downstream_b = focal_end+max_dist
            # if upstream_b <0:
                # upstream_b = 0
            # if downstream_b>len(mat):
                # downstream_b = len(mat)
            # upstream_region = mat[upstream_b:focal_start,upstream_b:focal_start]
            # downstream_region = mat[focal_end+1:downstream_b+1, focal_end+1:downstream_b+1]
            # mean_u = np.nanmean([np.exp(i) for i in list(chain(*upstream_region))])
            # mean_d = np.nanmean([np.exp(i) for i in list(chain(*upstream_region))])
            # insulation_track.append(max(mean_u, mean_d)/center_mean)
    # return(insulation_track)

def comparePreds(indiv1, indiv2):
    mse = np.mean(np.square(indiv1 - indiv2))
    spearman = stats.spearmanr(indiv1, indiv2)[0]

    # indiv1_triangle = [np.mean(indiv1[i]) for i in listOfIndexes]
    # indiv2_triangle = [np.mean(indiv2[i]) for i in listOfIndexes]
    # triangle_spearman = stats.spearmanr(indiv1_triangle, indiv2_triangle)[0]
    # triangle_mse = np.mean(np.square(np.array(indiv1_triangle) - np.array(indiv2_triangle)))
    # mat1 = from_upper_triu(indiv1, 448, 2)
    # mat2 = from_upper_triu(indiv2, 448, 2)
    # insulation_spearman = stats.spearmanr(get_insulation_track(mat1,10,100),get_insulation_track(mat2,10,100),nan_policy="omit")[0]

    # return (mse, spearman, triangle_mse, triangle_spearman, insulation_spearman)
    return (mse, spearman)

#####################################

# Useful for slurm scripts
indivname1=sys.argv[1].split(" ")[0]
indivname2=sys.argv[1].split(" ")[1]

#indivname1=sys.argv[1]
#indivname2=sys.argv[2]
#path2indiv1=sys.argv[3]
#path2indiv2=sys.argv[4]

import gzip

print("Indiv1 = %s, Indiv2 = %s" % (indivname1, indivname2),flush=True)


in_file_loc1 = '/wynton/group/capra/projects/modern_human_3Dgenome/data/akitaPreds/3dpreds/'
in_file_loc2 = '/wynton/group/capra/projects/modern_human_3Dgenome/data/akitaPreds/3dpreds/'

# if using gzip file read binary and use .decode() for l1 & l2 below (4x)
if path.exists("%s3dpreds_%s.txt.gz" % (in_file_loc1, indivname1)):
    f1 = gzip.open("%s3dpreds_%s.txt.gz" % (in_file_loc1, indivname1),"rb")
else:
    f1 = open("%s3dpreds_%s.txt" % (in_file_loc1, indivname1),"rb")

if path.exists("%s3dpreds_%s.txt.gz" % (in_file_loc2, indivname2)):
    f2 = gzip.open("%s3dpreds_%s.txt.gz" % (in_file_loc2, indivname2),"rb")
else:
    f2 = open("%s3dpreds_%s.txt" % (in_file_loc2, indivname2),"rb")

f_out = open("/wynton/group/capra/projects/modern_human_3Dgenome/data/pairwise/3dcomp_%s_vs_%s.txt" % (indivname1,indivname2),"w")

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
