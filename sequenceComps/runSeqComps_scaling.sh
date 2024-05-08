#!/bin/bash
#$ -N runSeqComps
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=24:00:00
#$ -l mem_free=10G
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runSeqComps/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runSeqComps/
#$ -t 1-3


echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
# source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d



COMP_LIST=/wynton/group/capra/projects/modern_human_3Dgenome/data/reference/lists/rep_indivs.txt

indivs=$(awk -v var="$SGE_TASK_ID" 'NR==var' ${COMP_LIST})
echo $indivs

PYSCRIPT=runSeqComparisons_scaling.py

python3 ${PYSCRIPT} --indivs "$indivs" --window_size_exp "$EXP" > /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runSeqComps/runSeqComps_"$SGE_TASK_ID".python.out
