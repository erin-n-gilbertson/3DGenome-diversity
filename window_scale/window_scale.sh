#!/bin/bash
#$ -N window_scale
#$ -M erin.gilbertson@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/window_scale.out
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/window_scale.err
#$ -l h_rt=72:00:00
#$ -l mem_free=20G
#$ -t 1-130


# load conda environment
echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d

echo "SGE_TASK_ID:  ${SGE_TASK_ID}"

#cd /wynton/group/capra/projects/modern_human_3Dgenome/data/reference/lists

indiv=$(awk -v var="$SGE_TASK_ID" 'NR==var' /wynton/group/capra/projects/modern_human_3Dgenome/data/reference/lists/rep_indivs.txt)
echo "indiv: ${indiv}"
echo "window exponent: ${EXP}"

python3 /wynton/group/capra/projects/modern_human_3Dgenome/bin/window_scale/window_scale.py --indiv "$indiv" --window_size_exponent "$EXP"