#!/bin/bash
#$ -N in_silico_mutagenesis_pts_specific_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/in_silico_mutagenesis_pts_specific_variants.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/in_silico_mutagenesis_pts_specific_variants.err
#$ -l h_rt=24:00:00
#$ -l mem_free=20G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basenji

# change directories
cd ../../data/in_silico_mutagenesis

# run
python3 ../../scripts/3_in_silico_mutagenesis/in_silico_mutagenesis.py --fasta ../reference_fasta/filtered_panTro6.fa --input ../lineage_specific_variants/pts_specific_variants_with_window.txt --start 0 --end 610 --out pts_variants.txt
