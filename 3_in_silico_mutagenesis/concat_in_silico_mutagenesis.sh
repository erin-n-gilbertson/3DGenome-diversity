#!/bin/bash
#$ -N concat_in_silico_mutagenesis
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/concat_in_silico_mutagenesis.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/concat_in_silico_mutagenesis.err
#$ -l h_rt=1:00:00
#$ -l mem_free=10G

# change directories
cd ../../data/in_silico_mutagenesis

# ppn clustering
head -n 1 ppn_pt_clustering_variants_1.txt > ppn_pt_clustering_window_variants.txt
for i in {1..930}; do tail -n +2 ppn_pt_clustering_variants_$i.txt >> ppn.tmp; done
sort -k1,1 -k2n,2 ppn.tmp >> ppn_pt_clustering_window_variants.txt
rm ppn.tmp

# pte
#head -n 1 pte_variants_1.txt > pte_specific_variants.txt
#for i in {1..130}; do tail -n +2 pte_variants_$i.txt >> pte.tmp; done
#sort -k1,1 -k2n,2 pte.tmp >> pte_specific_variants.txt
#rm pte.tmp

# ptv
#head -n 1 ptv_variants_1.txt > ptv_specific_variants.txt
#for i in {1..46}; do tail -n +2 ptv_variants_$i.txt >> ptv.tmp; done
#sort -k1,1 -k2n,2 ptv.tmp >> ptv_specific_variants.txt
#rm ptv.tmp
