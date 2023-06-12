#!/bin/bash
#$ -N retrieve_ptt_specific_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/retrieve_ptt_specific_variants.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/retrieve_ptt_specific_variants.err
#$ -l h_rt=24:00:00
#$ -l mem_free=20G

# change directories
cd ../../data/lineage_specific_variants

# set paths
script_path="../../scripts/3_in_silico_mutagenesis/retrieve_private_variants.py"
genotypes_path="../vcfs/one_zero_genotypes.txt"

# run
python3 "$script_path" --genotypes "$genotypes_path" --ids 2 8 11 14 21 24 29 37 39 40 43 46 49 52 53 56 --out ptt_specific_variants.txt
