#!/bin/bash
#$ -N diff_shuffles
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=55:00:00
#$ -l mem_free=10G
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/diff_shuffles/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/diff_shuffles/

#!/usr/bin/env python

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
source /wynton/home/capra/egilbertson/envs/akita/bin/activate

#python /wynton/group/capra/projects/modern_human_3Dgenome/bin/diffShuffles/differentiation_shuffles_all.py
python /wynton/group/capra/projects/modern_human_3Dgenome/bin/diffShuffles/differentiation_shuffles_nonAFR.py