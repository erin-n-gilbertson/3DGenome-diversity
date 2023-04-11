#!/bin/bash
#$ -N runExpectedDivergence
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=100:00:00
#$ -l mem_free=10G
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/expectedDiv/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/expectedDiv/

#!/usr/bin/env python

cd /wynton/group/capra/projects/modern_human_3Dgenome/bin/expectedDivergence

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d

MOD_INDIV="AMR_CLM_female_HG01119"

ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

CHR=${ARRY[$SGE_TASK_ID-1]}
echo "chr: ${CHR}"

python3 expected3DgenomeDivergence_oneChromOneInd.py "$MOD_INDIV" "$CHR"
