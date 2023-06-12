#!/bin/bash
#$ -N retrieve_pte_specific_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/retrieve_pte_specific_variants.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/retrieve_pte_specific_variants.err
#$ -l h_rt=24:00:00
#$ -l mem_free=20G

# change directories
cd ../../data/divergent_windows

# set paths
script_path="../../scripts/3_in_silico_mutagenesis/retrieve_private_variants.py"
genotypes_path="../vcfs/one_zero_genotypes.txt"

# run
python3 "$script_path" --genotypes "$genotypes_path" --ids 1 19 30 35 48 --out pte.tmp
