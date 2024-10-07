#!/bin/bash
#$ -N runAkitaIndiv
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=1:00:00
#$ -l mem_free=10G
#$ -t 1-4873


echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc

load_conda
conda activate modern3d
cd /wynton/group/capra/projects/modern_human_3Dgenome/EXAMPLE

chrm=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' ../windows.txt)
pos=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' ../windows.txt)

python example_individual_predictions.py --chrm $chrm --window_pos $pos 
