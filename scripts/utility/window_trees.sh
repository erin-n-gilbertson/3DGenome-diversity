#!/bin/bash
#$ -N window_trees
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/other_scripts/window_trees.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/other_scripts/window_trees.err
#$ -l h_rt=12:00:00
#$ -l mem_free=20G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate jupyter

# change directories
cd ../../data

# run
python3 ../scripts/other_scripts/window_trees.py
