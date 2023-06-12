#!/bin/bash
#$ -N IDW_3d_modifying_variant_prediction_and_comparison
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/IDW_3d_modifying_variant_prediction_and_comparison.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/IDW_3d_modifying_variant_prediction_and_comparison.err
#$ -l h_rt=2:00:00
#$ -l mem_free=30G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basenji

# change directories
cd ../../data/IDWs

# run
python ../../scripts/3_in_silico_mutagenesis/IDW_3d_modifying_variant_prediction_and_comparison.py --fasta ../reference_fasta/filtered_panTro6.fa --input IDW_3d_modifying_variants.txt
