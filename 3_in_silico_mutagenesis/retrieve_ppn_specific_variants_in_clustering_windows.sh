#!/bin/bash
#$ -N retrieve_ppn_specific_variants_in_clustering_windows
#$ -t 1-339
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/retrieve_ppn_specific_variants_in_clustering_windows.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/retrieve_ppn_specific_variants_in_clustering_windows.err
#$ -l h_rt=12:00:00
#$ -l mem_free=20G

# change directories
cd ../../data/clustering_windows

# set paths
script_path="../../scripts/3_in_silico_mutagenesis/retrieve_private_variants.py"
genotypes_path="../vcfs/one_zero_genotypes.txt"

# set variables
chr=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' ../clustering_windows/ppn_pt_clustering_windows.bed)
start=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' ../clustering_windows/ppn_pt_clustering_windows.bed)
end=$(awk -v row=$SGE_TASK_ID 'NR == row {print $3}' ../clustering_windows/ppn_pt_clustering_windows.bed)

# run
python3 "$script_path" --genotypes "$genotypes_path" --ids 9 20 22 25 26 33 34 36 45 --out "$chr"_"$start".tmp --region "$chr" $(($start+1)) $(($end+1))
