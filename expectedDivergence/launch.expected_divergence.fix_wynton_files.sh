#!/bin/bash
#$ -N runExpectedDivergence
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=100:00:00
#$ -l mem_free=10G
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/expectedDiv_fix/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/expectedDiv_fix/

#!/usr/bin/env python

cd /wynton/group/capra/projects/modern_human_3Dgenome/bin/expectedDivergence

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d

INDIVS=("AFR_ESN_female_HG03105" "AFR_ESN_female_HG03105" "AMR_CLM_female_HG01119" "AMR_CLM_female_HG01119" "EAS_CDX_female_HG00759" "SAS_BEB_female_HG03007" "SAS_BEB_female_HG03007")

ARRY=(16 19 19 5 15 2 7)

CHR=${ARRY[$SGE_TASK_ID-1]}
MOD_INDIV=${INDIVS[$SGE_TASK_ID-1]}
echo "ind: ${MOD_INDIV}"
echo "chr: ${CHR}"

python3 expected3DgenomeDivergence_oneChromOneInd.py "$MOD_INDIV" "$CHR"
