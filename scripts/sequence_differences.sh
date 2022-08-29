#!/bin/bash
#$ -N sequence_differences
#$ -t 1-2485
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/sequence_differences.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/sequence_differences.err
#$ -l h_rt=24:00:00
#$ -l mem_free=20G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ancestral_allele

# change directories
cd ../data/sequence_differences

# assign variables using the SGE task ID to get a specific pair of individuals
ind1=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' ../sample_pairs.txt)
ind2=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' ../sample_pairs.txt)

# run
python ../../scripts/sequence_differences.py --intervals ../pantro6_intervals.txt --sample_1 ../fastas/"$ind1".fa --sample_2 ../fastas/"$ind2".fa --sample_1_id "$ind1" --sample_2_id "$ind2"
