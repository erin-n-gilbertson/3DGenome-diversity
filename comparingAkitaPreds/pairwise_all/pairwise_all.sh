#!/bin/bash
#$ -N pairwise_all
#$ -M erin.gilbertson@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/
#$ -l h_rt=144:00:00
#$ -l mem_free=100G
#$ -t 1-22


# load conda environment
echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d


cd /wynton/group/capra/projects/modern_human_3Dgenome
pwd
ARRY=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)
CHR=${ARRY[$SGE_TASK_ID-1]}

echo "call python"
python3 bin/pairwise_all/pairwise_all.py --chromosome "$CHR" > stdout/pairwise_all_"$CHR".python.out