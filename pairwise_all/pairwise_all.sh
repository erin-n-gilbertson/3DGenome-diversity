#!/bin/bash
#$ -N pairwise_all
#$ -M erin.gilbertson@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/pairwise_all.out
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/pairwise_all.err
#$ -l h_rt=144:00:00
#$ -l mem_free=100G


# load conda environment
echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d


cd /wynton/group/capra/projects/modern_human_3Dgenome

echo "call python"
python3 bin/pairwise_all/pairwise_all.py > /wynton/group/capra/projects/modern_human_3Dgenome/stdout/pairwise_all.python.out